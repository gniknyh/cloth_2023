#include <mitsuba/core/fwd.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/string.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/microfacet.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/sampler.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class Sadeghi final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture, MicrofacetDistribution)


    Sadeghi(const Properties& props) : Base(props) {
        m_reflectance = props.texture<Texture>("reflectance", .5f);
        m_A = props.get<float>("albedo_volume", 0.3);
        m_eta = props.get<float>("eta", 1.5);
        m_inv_eta = 1.0f / m_eta;
        m_kd = props.get<float>("kd", 0.3);
        m_gamma_s = dr::deg_to_rad(props.get<float>("gamma_s", 12));
        m_gamma_v = dr::deg_to_rad(props.get<float>("gamma_v", 24));
        m_alpha = props.get<float>("alpha", 0.33);

 
        m_uv_scale = props.get<float>("uv_scale", 1.f);
        m_sample_reflect = true;
        m_sample_vol = true;

        BSDFFlags extra = BSDFFlags(0);
        m_components.push_back(BSDFFlags::GlossyReflection | BSDFFlags::FrontSide |
            BSDFFlags::BackSide | extra);
        m_components.push_back(BSDFFlags::GlossyTransmission | BSDFFlags::FrontSide |
            BSDFFlags::BackSide | BSDFFlags::NonSymmetric | extra);
        m_flags = m_components[0] | m_components[1];
        dr::set_attr(this, "flags", m_flags);

        parameters_changed();
    }

    void traverse(TraversalCallback* callback) override {
        callback->put_parameter("eta", m_eta, ParamFlags::Differentiable | ParamFlags::Discontinuous);
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/ = {}) override {
    } 

    // coord 
    /* extract theta coordinate from 3D direction
     * -pi < theta < pi */
    Float dir_theta(const Vector3f& w) const
    {
        return dr::asin(w.x());
    }

    /* extract phi coordinate from 3D direction.
     * -pi < phi < pi
     * Assuming phi(wi) = 0 */
    Float dir_phi(const Vector3f& w) const
    { 
        return dr::atan2(w.z(), w.y()); 
    }

    /* extract theta and phi coordinate from 3D direction
     * -pi/2 < theta < pi/2, -pi < phi < pi
     * Assuming phi(wi) = 0 */
    std::pair<Float, Float> dir_sph(const Vector3f& w) const
    {
        return std::make_pair(dir_theta(w), dir_phi(w));
    }

    /* compute the vector direction given spherical coordinates */ 
    Vector3f sph_dir(Float theta, Float phi) const
    {
        auto [sin_theta, cos_theta] = dr::sincos(theta);
        auto [sin_phi, cos_phi] = dr::sincos(phi);
        return Vector3f(sin_theta, cos_phi * cos_theta, cos_theta * sin_phi);
    }

    Float cauchy_pdf(Float gamma, Float x) const
    {
        Float value = dr::InvPi<Float> * gamma / (dr::sqr(x) + dr::sqr(gamma));
        return value;
    }

    Float normalized_gaussian(Float beta, Float theta) const{
		return dr::exp(-theta * theta / (2.0f * beta * beta)) / (dr::safe_sqrt(dr::TwoPi<Float> * beta * beta));
	}

    // t(0) = 1 and t(φ)monotonically decreases to almost zero as |φ|increases
    Float unit_height_gaussian(Float phi) const
    {
        // standard deviation chosen to be 20 deg
        Float sigma = dr::deg_to_rad(20.0f);

        Float denom = dr::safe_sqrt(2 * dr::Pi<Float>) * sigma;
        Float val = dr::exp(-dr::sqr(phi) / (2 * dr::sqr(sigma))) / denom;
        return val;
    }

    Float shadowing(Float phi_i, Float phi_o, Float phi_d) const
    {
        // Float  = phi_i - phi_o;
        Float u = unit_height_gaussian(phi_d);
        // u = dr::clamp(u, 0, 1);
            
        Float M_i = dr::maximum(0.f, dr::cos(phi_i));
        Float M_o = dr::maximum(0.f, dr::cos(phi_o));

        Float correlated = u * dr::minimum(M_i, M_o);
        Float uncorrelated = (1-u) * M_i * M_o;
        // M term
        Float M = correlated + uncorrelated;
        return M;
    }

    SurfaceInteraction3f scaled_uv_si(const SurfaceInteraction3f& si) const
    {
        SurfaceInteraction3f tmp(si);
        tmp.uv *= m_uv_scale;
        return tmp;
    }

    std::pair<Float, Float> sample_reflect(const Vector3f& w, const Point2f& sample2) const
    {
        Float xi_1 = sample2.x();
        Float xi_2 = sample2.y();

        auto [theta_r, phi_r] = dir_sph(w);
        auto A = dr::atan((theta_r + dr::Pi<Float> / 2) / (2 * m_gamma_s));
        auto B = dr::atan((theta_r - dr::Pi<Float> / 2) / (2 * m_gamma_s));
        Float theta_h = dr::tan(xi_1 * (A - B) + B);

        Float theta_i = 2 * m_gamma_s * theta_h - theta_r;
        Float phi_d = 2 * dr::asin(2 * xi_2 - 1);
        Float phi_i = phi_r + phi_d;
        
        return {theta_i, phi_i}; 
    }

    // std::pair<Float, Float> sample_volume_scatter(const Vector3f &w, const Point2f& sample2)
    std::pair<Float, Float> sample_volume_scatter(const Vector3f& w, const Point2f& sample2) const
    {
        Float xi_1 = sample2.x();
        Float xi_2 = sample2.y();

        auto [theta_i, phi_i] = dir_sph(w);
        Float theta_r, phi_r;
        auto A = dr::atan((theta_i + dr::Pi<Float> / 2) / (2 * m_gamma_v));
        auto B = dr::atan((theta_i - dr::Pi<Float> / 2) / (2 * m_gamma_v));
        Float theta_h = dr::tan(xi_1 * (A - B) + B);
        Float theta_r_ani = 2 * m_gamma_v * theta_h - theta_i;

        Float theta_r_iso = dr::asin(2 * xi_2 - 1);
        // Float pdf_iso = dr::cos(theta_r_iso / 2);
        // ani
        if (dr::any_or<true>(xi_1 < m_kd))
            theta_r = theta_r_ani;
        // iso
        else
            theta_r = theta_r_iso;

 
        phi_r = dr::Pi<Float> + phi_i;

        return {theta_r, phi_r};

    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext& ctx,
        const SurfaceInteraction3f& si,
        Float sample1,
        const Point2f& sample2,
        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Vector3f wi = dr::normalize(si.wi); 
        auto [theta_i, phi_i] = dir_sph(wi);

        BSDFSample3f bs = dr::zeros<BSDFSample3f>(); 
        if (unlikely(dr::none_or<false>(active)))
            return { bs, 0.f };
        // Determine the type of interaction
        bool has_reflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 0),
            has_transmission = ctx.is_enabled(BSDFFlags::GlossyTransmission, 1);

        if (has_reflection | has_transmission)
            ;
 
        Mask selected_r, selected_t;

        Float prob_reflectance = 0.5;
        // Float prob_reflectance = 1;
        selected_r = sample1 <= prob_reflectance && active;
        selected_t = !selected_r && active; 

        Float f;
        if (dr::any_or<true>(selected_r)) {
            bs.sampled_component = 0;
            bs.sampled_type = +BSDFFlags::GlossyReflection;
            Float theta_o, phi_o;
            std::tie(theta_o, phi_o) = sample_reflect(wi, sample2);
            Vector3f wo = sph_dir(theta_o, phi_o);
            bs.wo = wo; 
            bs.pdf = pdf_reflectance(wi, bs.wo );
            f = eval_reflectance(ctx, si, bs.wo, active); 
        }

        if (dr::any_or<true>(selected_t)) {
            bs.sampled_component = 1;
            bs.sampled_type = +BSDFFlags::GlossyTransmission;
            Float theta_o, phi_o;
            std::tie(theta_o, phi_o) = sample_volume_scatter(wi, sample2); 
            Vector3f wo = sph_dir(theta_o, phi_o);
            bs.wo = wo;
            bs.pdf = pdf_volume(wi, bs.wo );
            f = eval_volume(ctx, si, bs.wo, active);
        }

        

        auto [theta_o, phi_o] = dir_sph(bs.wo);
        active &= bs.pdf > 0.f;
        Spectrum r = m_reflectance->eval(scaled_uv_si(si), active);

        Float cos_theta_i = dr::cos(theta_i);
        cos_theta_i = dr::select(cos_theta_i > 0.f, cos_theta_i, 0.f);
        UnpolarizedSpectrum result(f * r * cos_theta_i);

        return { bs, (depolarizer<Spectrum>(result) / bs.pdf) & active };
    }

    Float eval_reflectance(const BSDFContext& ctx, const SurfaceInteraction3f& si,
        const Vector3f& wo, Mask active) const {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        bool has_reflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
        if (has_reflection)
            ;
 
        // fresnel
        Vector3f wi = dr::normalize(si.wi);
        Vector3f wo_norm = dr::normalize(wo);
        auto [theta_i, phi_i] = dir_sph(wi);
        auto [theta_o, phi_o] = dir_sph(wo);

        Vector3f wi_on_normal_plane = Vector3f(0, wi.y(), wi.z());
        Vector3f wo_on_normal_plane = Vector3f(0, wo_norm.y(), wo_norm.z());
        Float cos_phi_d = dr::dot(wi_on_normal_plane, wo_on_normal_plane);
        Float phi_d = dr::acos(cos_phi_d); 

        Float theta_d = (theta_i - theta_o) * 0.5f;
        Float theta_h = (theta_i + theta_o) * 0.5f;

        Float cos_theta_d = dr::cos(theta_d);
        Float cos_theta_d_2 = dr::sqr(cos_theta_d);
        Float cos_phi_d_half = dr::cos(phi_d / 2);

        Float F = 0.f;
        F = std::get<0>(fresnel(dr::cos(theta_d) * cos_phi_d_half, m_eta));

        // cauchy
        Float cauchy_term = cauchy_pdf(m_gamma_s, theta_h);
        // Float cauchy_term = normalized_gaussian(m_gamma_s, theta_h);
        // Float gau = normalized_gaussian(m_gamma_s, theta_h);
        
        Float x(0);
        x = F * cos_phi_d_half * cauchy_term;
        x /= cos_theta_d_2;

        if (dr::any_or<true>(cos_theta_d < 1e-5f))
            x = dr::clamp(x, 0.f, 1.f); 

        return x; 
    }

    Float eval_volume(const BSDFContext& ctx, const SurfaceInteraction3f& si,
        const Vector3f& wo, Mask active) const {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        // bool has_reflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
        // if (has_reflection)
        //     ; 

        // fresnel
        Vector3f wi = dr::normalize(si.wi);
        Vector3f wo_norm = dr::normalize(wo);
        auto [theta_i, phi_i] = dir_sph(wi);
        auto [theta_o, phi_o] = dir_sph(wo_norm);
        Float theta_h = (theta_i + theta_o) / 2;
        Float theta_d = (theta_i - theta_o) * 0.5f;

        Float cos_theta_d = dr::cos(theta_d);
        Float cos_theta_d_2 = dr::sqr(cos_theta_d);

        // fresnel
        Float cos_theta_i = dr::cos(theta_i);
        Float cos_theta_o = dr::cos(theta_o);

        auto [F_r_1, cos_theta_t, eta_it, eta_ti] =  fresnel(cos_theta_i, m_eta);
        Float F_t_1 = 1.0f - F_r_1;
        Float F_r_x = std::get<0>(fresnel(cos_theta_t, m_inv_eta));
        Float F_t_2 = 1.0f - F_r_x;
        Float F = F_t_1 * F_t_2;

        // cauchy
        Float cauchy_term = cauchy_pdf(m_gamma_v, theta_h);

        // cos
        Float denom = dr::cos(theta_i) + dr::cos(theta_o);
        active &= dr::neq(denom, 0.f);

        Float x =  F * m_A * ((1.0f - m_kd) * cauchy_term + m_kd) / denom;
        x /= cos_theta_d_2;

        if (dr::any_or<true>(cos_theta_d < 1e-5f))
            x = dr::clamp(x, 0.f, 1.f); 

        return x;
    }

    Spectrum eval(const BSDFContext& ctx, const SurfaceInteraction3f& si,
        const Vector3f& wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Spectrum result(0.f);
        bool has_reflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
        bool has_transmission = ctx.is_enabled(BSDFFlags::GlossyTransmission, 1);


        Vector3f wi = dr::normalize(si.wi);
        Vector3f wo_norm = dr::normalize(dr::normalize(wo));

        Float theta_i, phi_i;
        std::tie(theta_i, phi_i) = dir_sph(wi);
        auto [theta_o, phi_o] = dir_sph(wo_norm);

        Vector3f wi_on_normal_plane = Vector3f(0, wi.y(), wi.z());
        Vector3f wo_on_normal_plane = Vector3f(0, wo_norm.y(), wo_norm.z());
        Float cos_phi_d = dr::dot(wi_on_normal_plane, wo_on_normal_plane);
        Float phi_d = dr::acos(cos_phi_d);
        
        if (unlikely(dr::none_or<false>(active)))
            return result;

        Float M = shadowing(phi_i, phi_o, phi_d); 
        Float fr = 0, fv = 0;
        
        if (m_sample_reflect)
            fr += eval_reflectance(ctx, si, wo, active);
        if (m_sample_vol)
            fv += eval_volume(ctx, si, wo, active);

        Spectrum r = m_reflectance->eval(scaled_uv_si(si), active);
        Float f = 0.f;

        // maybe some bugs with shadowing
        // f = (fr+fv) * M ;
        f = (fr+fv);
        result = r * f; 
 
        return result & active;
    }

    Float pdf_reflectance(const Vector3f& wi, const Vector3f& wo) const
    {
        auto [theta_i, phi_i] = dir_sph(wi);
        auto [theta_o, phi_o] = dir_sph(wo);
        auto A = dr::atan((theta_i + dr::Pi<Float> / 2) / (2 * m_gamma_s));
        auto B = dr::atan((theta_i - dr::Pi<Float> / 2) / (2 * m_gamma_s));

        Float theta_h = (theta_i + theta_o) / 2;
        Float phi_d = phi_i - phi_o;
        phi_d = dr::select(phi_d > dr::Pi<Float>, phi_d - dr::Pi<Float>, phi_d);
        phi_d = dr::select(phi_d < -dr::Pi<Float>, phi_d + dr::Pi<Float>, phi_d);

        Float pdf_theta = 1.0f / (2.0f * (A - B) * dr::cos(theta_i)) * (m_gamma_s / (dr::sqr(theta_h) + dr::sqr(m_gamma_s)));
        Float pdf_phi = 0.25f * dr::cos(phi_d / 2.0f);
        Float pdf = pdf_theta * pdf_phi;
        return pdf;
    }

    Float pdf_volume(const Vector3f& wi, const Vector3f& wo) const
    {
        auto [theta_i, phi_i] = dir_sph(wi);
        auto [theta_o, phi_o] = dir_sph(wo);

        Float theta_h = (theta_i + theta_o) / 2;

        auto A = dr::atan((theta_i + dr::Pi<Float> / 2) / (2 * m_gamma_v));
        auto B = dr::atan((theta_i - dr::Pi<Float> / 2) / (2 * m_gamma_v));

        Float pdf_ani = 1.0f / (2.0f * (A - B) * dr::cos(theta_i)) * (m_gamma_v / (dr::sqr(theta_h) + dr::sqr(m_gamma_v)));

        Float pdf_iso = dr::cos(theta_o) / 2;
        Float pdf = (1 - m_kd) * pdf_ani + m_kd * pdf_iso;

        return pdf;
    }


    Float pdf(const BSDFContext& ctx, const SurfaceInteraction3f& si,
        const Vector3f& wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        // bool has_reflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
        // if (has_reflection)
        //     ;
        Vector3f wi = dr::normalize(si.wi);
        Float pdf_r = pdf_reflectance(wi, wo);
        Float pdf_v = pdf_volume(wi, wo);
        Float fresnel_p = std::get<0>(fresnel(Frame3f::cos_theta(wi), m_eta));

        Float p = 0;
        if (m_sample_reflect && m_sample_vol)
            p =  pdf_r * fresnel_p + pdf_v * (1-fresnel_p);
        else if (m_sample_reflect)
            p = pdf_r ;
        else if (m_sample_vol)
            p = pdf_v ;


        return dr::select(active, p, 0.f);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "  eta = " << m_eta << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:


    ref<Texture> m_reflectance;
    Float m_A;   // color albedo coeffi
    Float m_eta, m_inv_eta;
    Float m_kd;
    Float m_gamma_s, m_gamma_v;
    Float m_alpha;

    Float m_uv_scale;
    bool m_sample_reflect;
    bool m_sample_vol;

};

MI_IMPLEMENT_CLASS_VARIANT(Sadeghi, BSDF)
MI_EXPORT_PLUGIN(Sadeghi, "cloth")
NAMESPACE_END(mitsuba)

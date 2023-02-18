#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/string.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. tabs::
    .. code-tab:: xml
        :name: ClothMap

        <bsdf type="ClothMap">
            <texture name="ClothMap" type="bitmap">
                <boolean name="raw" value="true"/>
                <string name="filename" value="textures/ClothMap.jpg"/>
            </texture>
            <bsdf type="roughplastic"/>
        </bsdf>
*/

template <typename Float, typename Spectrum>
class ClothMap final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    ClothMap(const Properties &props) : Base(props) {
        for (auto &[name, obj] : props.objects(false)) {
            auto bsdf = dynamic_cast<Base *>(obj.get());

            if (bsdf) {
                if (m_nested_bsdf)
                    Throw("Only a single BSDF child object can be specified.");
                m_nested_bsdf = bsdf;
                props.mark_queried(name);
            }
        }
        if (!m_nested_bsdf)
            Throw("Exactly one BSDF child object must be specified.");

        // TODO: How to assert this is actually a RGBDataTexture?
        m_uv_scale = props.get<float>("uv_scale", 1.0f);
        m_normalmap = props.texture<Texture>("normalmap");
        m_tangentmap = props.texture<Texture>("tangentmap");

        // Add all nested components
        m_flags = (uint32_t) 0;
        for (size_t i = 0; i < m_nested_bsdf->component_count(); ++i) {
            m_components.push_back((m_nested_bsdf->flags(i)));
            m_flags |= m_components.back();
        }
        dr::set_attr(this, "flags", m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("nested_bsdf", m_nested_bsdf.get(), +ParamFlags::Differentiable);
        callback->put_object("NormalMap",   m_normalmap.get(),   ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("TangentMap",   m_tangentmap.get(),   ParamFlags::Differentiable | ParamFlags::Discontinuous);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float sample1,
                                             const Point2f &sample2,
                                             Mask active) const override {
        // Sample nested BSDF with perturbed shading frame
        SurfaceInteraction3f perturbed_si(si);
        perturbed_si.sh_frame = frame(si, active);
        perturbed_si.wi = perturbed_si.to_local(si.wi);
        auto [bs, weight] = m_nested_bsdf->sample(ctx, perturbed_si,
                                                  sample1, sample2, active);
 


        active &= dr::any(dr::neq(unpolarized_spectrum(weight), 0.f));
        if (dr::none_or<false>(active))
            return { bs, 0.f };

        // Transform sampled 'wo' back to original frame and check orientation
        Vector3f perturbed_wo = perturbed_si.to_world(bs.wo); 
        bs.wo = perturbed_wo;

        return { bs, weight & active };
    }


    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        // Evaluate nested BSDF with perturbed shading frame
        SurfaceInteraction3f perturbed_si(si);
        perturbed_si.sh_frame = frame(si, active);
        perturbed_si.wi       = perturbed_si.to_local(si.wi);
        Vector3f perturbed_wo = perturbed_si.to_local(wo);
 
        return m_nested_bsdf->eval(ctx, perturbed_si, perturbed_wo, active) & active;
    }


    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        // Evaluate nested BSDF with perturbed shading frame
        SurfaceInteraction3f perturbed_si(si);
        perturbed_si.sh_frame = frame(si, active);
        perturbed_si.wi       = perturbed_si.to_local(si.wi);
        Vector3f perturbed_wo = perturbed_si.to_local(wo);

        return dr::select(active, m_nested_bsdf->pdf(ctx, perturbed_si, perturbed_wo, active), 0.f);
    }
    
    SurfaceInteraction3f scaled_uv_si(const SurfaceInteraction3f& si) const
    {
        SurfaceInteraction3f tmp(si);
        tmp.uv *= m_uv_scale;
        return tmp;
    }

    Frame3f frame(const SurfaceInteraction3f &si, Mask active) const {
        SurfaceInteraction3f si_scaled = scaled_uv_si(si);
        Normal3f n = dr::fmadd(m_normalmap->eval_3(si_scaled, active), 2, -1.f);
        Vector3f t = dr::fmadd(m_tangentmap->eval_3(si_scaled, active), 2, -1.f);

        Frame3f result;
        result.n = dr::normalize(n);
        result.t = dr::normalize(t);
        result.s = dr::cross(result.n, result.t);
 
        return result;
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "ClothMap[" << std::endl
            << "  nested_bsdf = " << string::indent(m_nested_bsdf) << ","
            << std::endl
            << "  ClothMap = " << string::indent(m_normalmap) << ","
            << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
protected:
    ref<Base> m_nested_bsdf;
    ref<Texture> m_normalmap;
    ref<Texture> m_tangentmap;
    Float m_uv_scale;

};

MI_IMPLEMENT_CLASS_VARIANT(ClothMap, BSDF)
MI_EXPORT_PLUGIN(ClothMap, "Normal map material adapter");
NAMESPACE_END(mitsuba)

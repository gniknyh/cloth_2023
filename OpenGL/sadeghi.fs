#version 330 core
out vec4 FragColor;

struct Cloth
{
    float eta;
    float kd;
    float gamma_s;
    float gamma_v;
    float albedo;
};

float PI = 3.14159265359;

float lightRadiance = 25.0f;


in vec3 FragPos;  
in vec3 vertexNormal;  
in vec2 TexCoords;
in mat3 vs_TBN;

uniform sampler2D _tex;
uniform sampler2D normalMap;
uniform sampler2D tangentMap;

uniform vec3 viewPos;   // camera world position
uniform vec3 lightPos;  // light world position

uniform Cloth cloth;
 
// gaussian
float normalized_gaussian(float theta, float beta) {
    return exp(-theta * theta / (2.0f * beta * beta)) / (sqrt(2.0f * PI * beta * beta)); 
}
 
float dir_theta(vec3 w)
{
    return asin(w.x);
}

// float dir_phi(vec3 w)
// {
//     return atan(w.z, w.y); 
// }

float sqr(float x)
{
    return x * x;  
}


float fresnelSchlick(float cosTheta, float eta1, float eta2)
{
    float tmp = (eta1-eta2)/(eta1+eta2);
    float F0 = tmp*tmp;
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

float u_gaussian(float x) {
    // using 20 deg as std;
    float std = 0.349;
    float sqr_sd_2 = std * std * 2.0;
    return exp(-x * x / sqr_sd_2);
}

float shadowing(float cosPhiI, float cosPhiO, float phi_d) {
    float m_i = max(cosPhiI, 0.0);
    float m_o = max(cosPhiO, 0.0);
    float u = u_gaussian(phi_d);
    float corrated = min(m_i, m_o);
    float uncorrated = m_i * m_o;
    return mix(uncorrated, corrated, u);
}
  

float reflectance(float theta_i, float theta_r, float phi_d)
{

    // compute the longitudinal angles and azimuthal angle
    float theta_d = (theta_i - theta_r) * 0.5f;
    float theta_h = (theta_i + theta_r) * 0.5f;
    // float phi_d = phi_i - phi_r;
    float cos_theta_d_2 = cos(theta_d)*cos(theta_d);
    float cos_phi_d_half = cos(phi_d / 2);

    // cos_theta_d * cos_phi_d_half
    float fresnel_cos = cos(theta_d) * cos_phi_d_half;

    // fresnel term
    float fresnel_term = fresnelSchlick(fresnel_cos, cloth.eta, 1.f);
    // float fresnel_term = fresnelSchlick(cos(theta_i), cloth.eta, 1.f);

    // cosine term
    float cos_term = cos_phi_d_half;

    // gaussian term
    float gaussian_term = normalized_gaussian(theta_h, cloth.gamma_s);

    float result = 0.f;
    result = fresnel_term * cos_phi_d_half * gaussian_term;
    
    result /= cos_theta_d_2;

    return max(result, 0.f);
}

float volumetric_scatter(float theta_i, float theta_r, float phi_d)
{
    float theta_d = (theta_i - theta_r) * 0.5f;
    float theta_h = (theta_i + theta_r) * 0.5f;
    // float phi_d = phi_i - phi_r;
    float cos_theta_d_2 = cos(theta_d)*cos(theta_d);
    float cos_phi_d_half = cos(phi_d / 2);

    // fresnel term
    float f_r = fresnelSchlick(cos(theta_i), cloth.eta, 1.f);
    float fresnel_term = (1 - f_r) * (1 - f_r);

    // gaussian term
    float gaussian_term = normalized_gaussian(theta_h, cloth.gamma_v);

    float denom = cos(theta_i) + cos(theta_r);

    float result = fresnel_term * cloth.albedo * ((1-cloth.kd) * gaussian_term + cloth.kd) / denom;
    result /= cos_theta_d_2;

    return max(result, 0.f);
}

vec3 get_map_normal()
{
    vec3 normal = texture(normalMap, TexCoords).rgb;
    normal = normalize(normal * 2.0 - 1.0);  
    return normal;
}

vec3 get_map_tangent()
{
    vec3 tangent = texture(tangentMap, TexCoords).rgb;
    tangent = normalize(tangent * 2.0 - 1.0);  
    return tangent;
}

float brdf(vec3 wi, vec3 wo)
{
    wi = normalize(wi);
    wo = normalize(wo);

    // compute the longitudinal angles and azimuthal angle
    mat3 TBN = vs_TBN;
    vec3 normal = TBN * get_map_normal();
    vec3 tangent = TBN * get_map_tangent();

    // vec3 B = TBN * get_map_tangent();
    // vec3 tangent = cross(normal, B);
    
    float sin_theta_i = dot(tangent, wi);
    float theta_i = asin(sin_theta_i);

    float sin_theta_o = dot(tangent, wo);
    float theta_o = asin(sin_theta_o);

    vec3 wi_on_normal_plane = normalize(wi - tangent * sin_theta_i);
    vec3 wo_on_normal_plane = normalize(wo - tangent * sin_theta_o);
    float cos_phi_d = dot(wi_on_normal_plane, wo_on_normal_plane);
    float phi_d = acos(cos_phi_d);

    float cos_phi_i = dot(normal, wi_on_normal_plane);
    float cos_phi_o = dot(normal, wo_on_normal_plane);

    float shadow = shadowing(cos_phi_i, cos_phi_o, phi_d);
 
    // compute reflectance term
    float r = reflectance(theta_i, theta_o, phi_d);

    // compute vol scatter
    float v = volumetric_scatter(theta_i, theta_o, phi_d);

    float result = r + v;
    result *= shadow;

    return result;
}

void main()
{ 
    FragColor = vec4(0.f);

    vec3 texColor;
    texColor = texture(_tex, TexCoords).rgb;

    vec3 light_dir = normalize(lightPos - FragPos);
    vec3 view_dir = normalize(viewPos - FragPos);

    vec3 world_light_dir = normalize(lightPos - FragPos);
    vec3 world_view_dir = normalize(viewPos - FragPos);

    // directional light
    // world_light_dir = normalize(vec3(0.f,1.f,0.f));

    float light_cos = dot(world_light_dir, vertexNormal);
    float distance    = length(lightPos - FragPos);
    float attenuation = 1.0 / (distance * distance);
    float radiance     = lightRadiance * attenuation;   
    
    if (light_cos > 0)
    {
        vec3 wi =  world_light_dir;
        vec3 wo =  world_view_dir;

        float result = brdf(wi, wo);
        result *= light_cos; 

        FragColor = vec4(radiance*result*texColor, 1.0f); 
    } 
    
} 
#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;
layout (location = 3) in vec3 aTangent;

out vec3 FragPos;
out vec3 vertexNormal;
out vec2 TexCoords;
out mat3 vs_TBN; 

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform vec3 viewPos;   // camera world position
uniform vec3 lightPos;  // light world position

void main()
{
    FragPos = vec3(model * vec4(aPos, 1.0));
    mat3 normalMat = mat3(transpose(inverse(model)));
    vertexNormal = normalize(normalMat * aNormal);  
    vec3 vertexTangent = normalize(normalMat * aTangent);
    vec3 vertexBitangent = normalize(cross(vertexNormal, vertexTangent));
    TexCoords = aTexCoords * 200.f;

    vs_TBN = (mat3(vertexTangent, vertexBitangent, vertexNormal));

    gl_Position = projection * view * vec4(FragPos, 1.0);
}
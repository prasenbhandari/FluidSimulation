#version 430

in vec3 vertexPosition;
in vec2 vertexTexCoord;

uniform mat4 mvp;

struct Particle {
    vec2 position;
    vec2 predicted_position;
    vec2 velocity;
    vec2 force;
    float density;
    float near_density;
    float pressure;
    float padding;
};

layout(std430, binding = 0) buffer ParticleBuffer {
    Particle particles[];
};

out vec2 fragTexCoord;
out float fragDensity;

void main() {
    Particle p = particles[gl_InstanceID];
    float particleSize = 0.2;

    vec2 scaledVertex = vertexPosition.xz * particleSize;

    vec2 worldPosition = p.position + scaledVertex;
    
    gl_Position = mvp * vec4(worldPosition, 0.0, 1.0);

    fragTexCoord = vertexTexCoord;
    fragDensity = p.density;
}
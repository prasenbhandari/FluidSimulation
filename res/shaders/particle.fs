#version 430

in vec2 fragTexCoord;
in float fragSpeed;

const float maxVelocity = 8.0;

out vec4 fragColor;

void main() {
    vec2 center = vec2(0.5, 0.5);
    float dist = distance(fragTexCoord, center);
    
    float alpha = 1.0 - smoothstep(0.4, 0.45, dist);
    
    if (alpha <= 0.0) discard;
    
    float t = clamp(fragSpeed / maxVelocity, 0.0, 1.0);
    
    vec3 colorBlue = vec3(0.13, 0.34, 0.73);
    vec3 colorCyan = vec3(0.30, 1.0, 0.56);
    vec3 colorYellow = vec3(1.0, 0.93, 0.0);
    vec3 colorRed = vec3(0.97, 0.29, 0.03);
    
    vec3 color;
    if (t < 0.5) {
        color = mix(colorBlue, colorCyan, t * 2.0);
    } else if (t < 0.7) {
        color = mix(colorCyan, colorYellow, (t - 0.5) * 5.0);
    } else {
        color = mix(colorYellow, colorRed, (t - 0.7) * 3.33);
    }
    
    fragColor = vec4(color, alpha);
}

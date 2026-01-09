#version 430

in vec2 fragTexCoord;
in float fragDensity;

out vec4 fragColor;

void main() {
    vec2 center = vec2(0.5, 0.5);
    float dist = distance(fragTexCoord, center);
    
    float alpha = 1.0 - smoothstep(0.4, 0.45, dist);
    
    if (alpha <= 0.0) discard;
    
    float normalizedDensity = clamp(fragDensity / 20.0, 0.0, 1.0); 
    
    vec3 blue = vec3(0.0, 0.0, 1.0);
    vec3 cyan = vec3(0.0, 1.0, 1.0);
    vec3 white = vec3(1.0, 1.0, 1.0);
    
    vec3 color;
    if (normalizedDensity < 0.5) {
        color = mix(blue, cyan, normalizedDensity * 2.0);
    } else {
        color = mix(cyan, white, (normalizedDensity - 0.5) * 2.0);
    }
    
    fragColor = vec4(color, alpha);
}

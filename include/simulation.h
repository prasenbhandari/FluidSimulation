#pragma once

#include <vector>

#include "compute_shader.h"
#include "raylib.h"

const float PIXELS_PER_METER = 50.0f;

// Precomputed kernel constants
struct KernelConstants {
    float poly6_scale;            // 4 / (PI * h^8)
    float spiky_pow2_scale;       // 6 / (PI * h^4) for density
    float spiky_pow3_scale;       // 10 / (PI * h^5) for near density
    float spiky_pow2_grad_scale;  // -12 / (PI * h^4) for pressure gradient
    float spiky_pow3_grad_scale;  // -30 / (PI * h^5) for near pressure gradient
    float viscosity_scale;        // 40 / (PI * h^5)
    float smoothing_radius_sq;    // h^2
    float inv_smoothing_radius;   // 1 / h
};

struct SimParams {
    float gravity = 10.0f;
    float mass = 1.0f;
    float rest_density = 5.0f;
    float stiffness = 50.0f;
    float near_pressure_multiplier = 1.0f;
    float smoothing_radius = 0.5f;
    float damping = 0.97f;
    float boundary_damping = 0.95f;
    float interaction_radius = 2.0f;
    float interaction_strength = 20.0f;
    float viscosity_strength = 0.3f;
    bool gravity_enabled = true;
    int num_particles = 2000;
};

struct Particle {
    Vector2 position;
    Vector2 predicted_position;
    Vector2 velocity;
    Vector2 force;
    float density;
    float near_density;
    float pressure;
    float padding;
};

struct ParticleIndex {
    unsigned int cell_key;
    unsigned int particle_id;
};

class FluidSimulation {
public:
    FluidSimulation();
    ~FluidSimulation();

    void init();
    void reset();

    void update(float delta_time);
    void draw();
    void draw_gui();

    SimParams params;

private:
    unsigned int ssbo_indices;
    unsigned int ssbo_offsets;
    unsigned int ssbo_id;
    ComputeShader offset_shader;
    ComputeShader density_shader;
    ComputeShader force_shader;
    Shader particle_shader;
    ComputeShader integrate_shader;
    ComputeShader hash_shader;
    Mesh quad_mesh;
    Material quad_material;
    std::vector<Matrix> dummy_transforms;

    KernelConstants kernel_constants;
    float cached_smoothing_radius = 0.0f;

    void apply_mouse_force();
    void update_kernel_constants();

    void update_integration(float delta_time);
    void update_spatial_hash();
    void update_density();
    void update_forces();
};

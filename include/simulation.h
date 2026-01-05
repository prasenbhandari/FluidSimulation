#pragma once

#include <vector>

#include "raylib.h"
#include "spatial_hash.h"

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
    float rest_density = 10.0f;
    float stiffness = 50.0f;
    float near_pressure_multiplier = 2.0f;
    float smoothing_radius = 0.6f;
    float damping = 0.97f;
    float boundary_damping = -0.5f;
    float interaction_radius = 2.0f;
    float interaction_strength = 1200.0f;
    float viscosity_strength = 0.5f;
    bool gravity_enabled = true;
    int num_particles = 800;
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
    std::vector<Vector2> position;
    std::vector<Vector2> predicted_position;
    std::vector<Vector2> velocity;
    std::vector<Vector2> force;
    std::vector<float> density;
    std::vector<float> near_density;
    std::vector<float> pressure;

    SpatialHash spatial_hash;
    KernelConstants kernel_constants;
    float cached_smoothing_radius = 0.0f;

    void apply_mouse_force();
    void update_kernel_constants();

    void compute_density_pressure();
    void compute_forces();
    void compute_viscosity();
    void integrate(float delta_time);

    float poly6_kernel(float distance_squared);
    float spiky_pow2_kernel(float distance);
    float spiky_pow3_kernel(float distance);
    float spiky_pow2_gradient(float distance);
    float viscosity_kernel(float distance);
};

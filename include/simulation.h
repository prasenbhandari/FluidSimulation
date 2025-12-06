#pragma once

#include <vector>
#include "raylib.h"
#include "particle.h"
#include "spatial_hash.h"

const float PIXELS_PER_METER = 50.0f;

// Simulation parameters that can be adjusted at runtime
struct SimParams {
    float gravity = 10.0f;
    float mass = 1.0f;
    float rest_density = 10.0f;
    float stiffness = 50.0f;
    float smoothing_radius = 0.6f;
    float damping = 0.97f;
    float boundary_damping = -0.5f;
    float interaction_radius = 2.0f;
    float interaction_strength = 50.0f;
    bool gravity_enabled = false;
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
    void draw_gui();  // New method for raygui controls

    // Public parameters for GUI access
    SimParams params;

private:
    std::vector<Particle> particles;
    SpatialHash spatial_hash;

    void handle_input();

    // physics steps
    void compute_density_pressure();
    void compute_forces();
    void integrate(float delta_time);

    // SPH kernel functions
    float poly6_kernel(float distance_squared);
    float spiky_gradient(float distance);
};  

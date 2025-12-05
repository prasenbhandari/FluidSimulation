#pragma once

#include <vector>
#include "raylib.h"
#include "particle.h"

class FluidSimulation {
public:
    FluidSimulation();
    ~FluidSimulation();

    void init();
    void reset();

    void update(float delta_time);
    void draw();
    void handle_input();

private:
    std::vector<Particle> particles;

    // physics steps
    void compute_density_pressure();
    void compute_forces();
    void integrate(float delta_time);

    // SPH kernel functions
    float poly6_kernel(float distance_squared);
    float spiky_gradient(float distance);
};  

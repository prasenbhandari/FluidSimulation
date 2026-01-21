# Fluid Simulation

A fluid simulation using particle-based methods (Smoothed Particle Hydrodynamics) implemented in C++. The simulation runs on the GPU using compute shaders and is visualized with Raylib.

## Features

- **GPU-Accelerated Computation**: Uses OpenGL 4.3 compute shaders for fast particle simulations
- **Particle-Based Fluid Dynamics**: Implements SPH (Smoothed Particle Hydrodynamics) physics

## Requirements

- CMake 3.14 or higher
- C++17 compatible compiler
- OpenGL 4.3 capable GPU
- Raylib (automatically fetched by CMake)

## Building

```bash
mkdir build
cd build
cmake ..
make
```

## Running

```bash
./build/FluidSimRaylib
```

## Project Structure

- `src/` - Source files (main.cpp, simulation.cpp, compute_shader.cpp)
- `include/` - Header files (simulation.h, compute_shader.h, particle.h)
- `res/shaders/` - GLSL compute and rendering shaders

## Configuration

Key simulation parameters can be modified in `include/simulation.h`:

- `gravity` - Gravitational acceleration (default: 10.0f)
- `mass` - Particle mass (default: 1.0f)
- `rest_density` - Target fluid density (default: 5.0f)
- `smoothing_radius` - Kernel smoothing radius (default: 1.2f)
- `num_particles` - Number of particles to simulate (default: 400)

## How It Works

1. **Particle System**: Particles represent fluid elements with position, velocity, and density
2. **Density Calculation**: Compute shaders calculate particle densities using kernel functions
3. **Force Computation**: Pressure and viscosity forces are computed based on density
4. **Integration**: Velocities and positions are updated using physics integration
5. **Rendering**: Particles are rendered as circles at their positions

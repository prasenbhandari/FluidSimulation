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
./FluidSimRaylib
```

## Project Structure

- `src/` - Source files (main.cpp, simulation.cpp, compute_shader.cpp)
- `include/` - Header files (simulation.h, compute_shader.h, particle.h)
- `res/shaders/` - GLSL compute and rendering shaders

## References

- [Coding Adventure: Simulating Fluids](https://youtu.be/rSKMYc1CQHE?si=870lmoILmriTKVtd) by Sebastian Lague
- [Particle-Based Fluid Simulation for Interactive Applications](https://matthias-research.github.io/pages/publications/sca03.pdf) by Müller et al.
- [Smoothed Particle Hydrodynamics Techniques](https://sph-tutorial.physics-simulation.org/pdf/SPH_Tutorial.pdf) by Koschier et al.

#include "simulation.h"

#include <rlgl.h>

#include <algorithm>
#include <iostream>

#include "raymath.h"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

const int WINDOW_WIDTH = 800 / PIXELS_PER_METER;
const int WINDOW_HEIGHT = 600 / PIXELS_PER_METER;

FluidSimulation::FluidSimulation() {
    init();
    reset();
    update_kernel_constants();
}

FluidSimulation::~FluidSimulation() {}

void FluidSimulation::init() {
    ssbo_id = rlLoadShaderBuffer(sizeof(Particle) * params.num_particles, nullptr, RL_DYNAMIC_COPY);
    ssbo_indices = rlLoadShaderBuffer(sizeof(ParticleIndex) * params.num_particles, nullptr, RL_DYNAMIC_COPY);

    // Binding 2: Counts (Histogram)
    ssbo_counts = rlLoadShaderBuffer(sizeof(unsigned int) * params.num_particles, nullptr, RL_DYNAMIC_COPY);

    // Binding 3: Block Sums (Auxiliary for Scan)
    // Size = num_particles / 2048 (rounded up)
    int num_blocks = (params.num_particles + 2047) / 2048;
    ssbo_block_sums = rlLoadShaderBuffer(sizeof(unsigned int) * num_blocks, nullptr, RL_DYNAMIC_COPY);

    // Binding 4: Sorted Indices (Output of scatter)
    ssbo_sorted_indices = rlLoadShaderBuffer(sizeof(ParticleIndex) * params.num_particles, nullptr, RL_DYNAMIC_COPY);

    ssbo_offsets = rlLoadShaderBuffer(sizeof(unsigned int) * params.num_particles, nullptr, RL_DYNAMIC_COPY);

    particle_shader = LoadShader("../res/shaders/particle.vs", "../res/shaders/particle.fs");
    particle_shader.locs[SHADER_LOC_MATRIX_MVP] = GetShaderLocation(particle_shader, "mvp");

    hash_shader.load("../res/shaders/calculate_hash.comp");

    quad_mesh = GenMeshPlane(1.0f, 1.0f, 1, 1);
    quad_material = LoadMaterialDefault();
    quad_material.shader = particle_shader;

    integrate_shader.load("../res/shaders/integrate.comp");
    offset_shader.load("../res/shaders/calculate_offset.comp");
    density_shader.load("../res/shaders/calculate_density.comp");
    force_shader.load("../res/shaders/calculate_forces.comp");
    count_shader.load("../res/shaders/count_frequencies.comp");
    scan_shader.load("../res/shaders/scan_counts.comp");
    scatter_shader.load("../res/shaders/scatter.comp");
}

void FluidSimulation::reset() {
    std::vector<Particle> particles(params.num_particles);
    dummy_transforms.resize(params.num_particles);

    int particle_per_row = (int)sqrt(params.num_particles);
    float spacing = params.smoothing_radius * 0.9f;

    float startX = (WINDOW_WIDTH - particle_per_row * spacing) / 2.0f;
    float startY = WINDOW_HEIGHT / 4.0f;

    for (int i = 0; i < params.num_particles; i++) {
        dummy_transforms[i] = MatrixIdentity();

        float x = startX + (i % particle_per_row) * spacing;
        float y = startY + (i / particle_per_row) * spacing;

        x += (float)GetRandomValue(-10, 10) / 100.0f;

        y += (float)GetRandomValue(-10, 10) / 100.0f;

        particles[i].position = {x, y};
        particles[i].predicted_position = {x, y};
        particles[i].velocity = {0, 0};
        particles[i].force = {0, 0};
        particles[i].density = 0;
        particles[i].near_density = 0;
        particles[i].pressure = 0;
    }

    rlUpdateShaderBuffer(ssbo_id, particles.data(), sizeof(Particle) * params.num_particles, 0);
    update_kernel_constants();
}

void FluidSimulation::update(float delta_time) {
    update_integration(delta_time);
    update_spatial_hash();
    update_density();
    update_forces();
}

void FluidSimulation::update_integration(float delta_time) {
    unsigned int program = integrate_shader.get_program_id();
    rlEnableShader(program);

    int dt_loc = rlGetLocationUniform(program, "deltaTime");
    int grav_loc = rlGetLocationUniform(program, "gravity");
    int damp_loc = rlGetLocationUniform(program, "damping");
    int mass_loc = rlGetLocationUniform(program, "mass");
    int windowSize_loc = rlGetLocationUniform(program, "windowSize");

    rlSetUniform(dt_loc, &delta_time, SHADER_UNIFORM_FLOAT, 1);

    float current_gravity = params.gravity_enabled ? params.gravity : 0.0f;
    rlSetUniform(grav_loc, &current_gravity, SHADER_UNIFORM_FLOAT, 1);

    float damping_val = fabs(params.boundary_damping);
    rlSetUniform(damp_loc, &damping_val, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(mass_loc, &params.mass, SHADER_UNIFORM_FLOAT, 1);

    float size[2] = {(float)WINDOW_WIDTH, (float)WINDOW_HEIGHT};
    rlSetUniform(windowSize_loc, size, SHADER_UNIFORM_VEC2, 1);

    rlBindShaderBuffer(ssbo_id, 0);

    int num_groups = (params.num_particles + 63) / 64;
    integrate_shader.dispatch(num_groups, 1, 1);

    rlDisableShader();
}

void FluidSimulation::update_spatial_hash() {
    int num_groups = (params.num_particles + 63) / 64;

    unsigned int hash_program = hash_shader.get_program_id();
    rlEnableShader(hash_program);
    int space_loc = rlGetLocationUniform(hash_program, "spacing");
    int num_loc = rlGetLocationUniform(hash_program, "num_particles");

    unsigned int num_particles_uint = (unsigned int)params.num_particles;
    rlSetUniform(space_loc, &params.smoothing_radius, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(num_loc, &num_particles_uint, SHADER_UNIFORM_UINT, 1);

    rlBindShaderBuffer(ssbo_id, 0);
    rlBindShaderBuffer(ssbo_indices, 1);

    hash_shader.dispatch(num_groups, 1, 1);

    rlDisableShader();  // Disable hash shader

    // DEBUG: Check Hash output
    static int hash_debug_frame = 0;
    if (hash_debug_frame++ == 0) {
        std::vector<ParticleIndex> debug_indices(params.num_particles);
        rlReadShaderBuffer(ssbo_indices, debug_indices.data(), debug_indices.size() * sizeof(ParticleIndex), 0);
        std::cout << "DEBUG: Pre-Sort Keys: ";
        for (int i = 0; i < std::min((int)params.num_particles, 20); ++i) std::cout << debug_indices[i].cell_key << " ";
        std::cout << std::endl;
    }
    std::vector<unsigned int> zero_counts(params.num_particles, 0);
    rlUpdateShaderBuffer(ssbo_counts, zero_counts.data(), zero_counts.size() * sizeof(unsigned int), 0);

    unsigned int count_prog = count_shader.get_program_id();
    rlEnableShader(count_prog);
    rlSetUniform(rlGetLocationUniform(count_prog, "num_particles"), &num_particles_uint, SHADER_UNIFORM_UINT, 1);
    rlBindShaderBuffer(ssbo_indices, 1);
    rlBindShaderBuffer(ssbo_counts, 2);

    int count_groups = (params.num_particles + 255) / 256;
    count_shader.dispatch(count_groups, 1, 1);
    rlDisableShader();

    unsigned int scan_prog = scan_shader.get_program_id();
    rlEnableShader(scan_prog);
    rlSetUniform(rlGetLocationUniform(scan_prog, "num_particles"), &num_particles_uint, SHADER_UNIFORM_UINT, 1);
    rlBindShaderBuffer(ssbo_counts, 2);
    rlBindShaderBuffer(ssbo_block_sums, 3);

    int scan_groups = (params.num_particles + 2047) / 2048;
    scan_shader.dispatch(scan_groups, 1, 1);
    rlDisableShader();

    unsigned int scatter_prog = scatter_shader.get_program_id();
    rlEnableShader(scatter_prog);
    rlSetUniform(rlGetLocationUniform(scatter_prog, "num_particles"), &num_particles_uint, SHADER_UNIFORM_UINT, 1);
    rlBindShaderBuffer(ssbo_indices, 1);
    rlBindShaderBuffer(ssbo_counts, 2);
    rlBindShaderBuffer(ssbo_sorted_indices, 4);

    int scatter_groups = (params.num_particles + 255) / 256;
    scatter_shader.dispatch(scatter_groups, 1, 1);
    rlDisableShader();

    std::vector<ParticleIndex> sorted_cpu(params.num_particles);
    rlReadShaderBuffer(ssbo_sorted_indices, sorted_cpu.data(), sorted_cpu.size() * sizeof(ParticleIndex), 0);
    rlUpdateShaderBuffer(ssbo_indices, sorted_cpu.data(), sorted_cpu.size() * sizeof(ParticleIndex), 0);

    // DEBUG: Print first 10 keys
    static int debug_frame = 0;
    if (debug_frame++ == 0) {
        std::cout << "DEBUG: Sorted Keys sample: ";
        for (int i = 0; i < std::min((int)params.num_particles, 20); ++i) {
            std::cout << sorted_cpu[i].cell_key << " ";
        }
        std::cout << std::endl;
    }

    std::vector<unsigned int> clear_offsets(params.num_particles, -1);
    rlUpdateShaderBuffer(ssbo_offsets, clear_offsets.data(), clear_offsets.size() * sizeof(unsigned int), 0);

    unsigned int offset_program = offset_shader.get_program_id();
    rlEnableShader(offset_program);

    int num_loc_off = rlGetLocationUniform(offset_program, "num_particles");
    rlSetUniform(num_loc_off, &num_particles_uint, SHADER_UNIFORM_UINT, 1);

    rlBindShaderBuffer(ssbo_indices, 0);
    rlBindShaderBuffer(ssbo_offsets, 1);

    offset_shader.dispatch(num_groups, 1, 1);

    rlDisableShader();
}

void FluidSimulation::update_density() {
    int num_groups = (params.num_particles + 63) / 64;

    unsigned int density_program = density_shader.get_program_id();
    rlEnableShader(density_program);

    int num_loc = rlGetLocationUniform(density_program, "num_particles");
    int poly6_loc = rlGetLocationUniform(density_program, "poly6_scale");
    int spiky_loc = rlGetLocationUniform(density_program, "spiky_pow3_scale");
    int smooth_loc = rlGetLocationUniform(density_program, "smoothing_radius_sq");
    int spacing_loc = rlGetLocationUniform(density_program, "spacing");
    int mass_loc = rlGetLocationUniform(density_program, "mass");

    unsigned int num_particles_uint = (unsigned int)params.num_particles;
    rlSetUniform(num_loc, &num_particles_uint, SHADER_UNIFORM_UINT, 1);
    rlSetUniform(poly6_loc, &kernel_constants.poly6_scale, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(spiky_loc, &kernel_constants.spiky_pow3_scale, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(smooth_loc, &kernel_constants.smoothing_radius_sq, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(spacing_loc, &params.smoothing_radius, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(mass_loc, &params.mass, SHADER_UNIFORM_FLOAT, 1);

    rlBindShaderBuffer(ssbo_id, 0);
    rlBindShaderBuffer(ssbo_indices, 1);
    rlBindShaderBuffer(ssbo_offsets, 2);

    density_shader.dispatch(num_groups, 1, 1);

    rlDisableShader();
}

void FluidSimulation::update_forces() {
    int num_groups = (params.num_particles + 63) / 64;

    unsigned int force_program = force_shader.get_program_id();
    rlEnableShader(force_program);

    int num_loc = rlGetLocationUniform(force_program, "num_particles");
    int spacing_loc = rlGetLocationUniform(force_program, "spacing");
    int smooth_loc = rlGetLocationUniform(force_program, "smoothing_radius");
    int stiffness_loc = rlGetLocationUniform(force_program, "stiffness");
    int near_pres_loc = rlGetLocationUniform(force_program, "near_pressure_multiplier");
    int visc_str_loc = rlGetLocationUniform(force_program, "viscosity_strength");
    int mass_loc = rlGetLocationUniform(force_program, "mass");

    int spiky_grad_loc = rlGetLocationUniform(force_program, "spiky_pow2_grad_scale");
    int near_grad_loc = rlGetLocationUniform(force_program, "spiky_pow3_grad_scale");
    int visc_loc = rlGetLocationUniform(force_program, "viscosity_scale");
    int rest_density_loc = rlGetLocationUniform(force_program, "rest_density");

    unsigned int num_particles_uint = (unsigned int)params.num_particles;
    rlSetUniform(num_loc, &num_particles_uint, SHADER_UNIFORM_UINT, 1);
    rlSetUniform(spacing_loc, &params.smoothing_radius, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(smooth_loc, &params.smoothing_radius, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(stiffness_loc, &params.stiffness, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(near_pres_loc, &params.near_pressure_multiplier, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(visc_str_loc, &params.viscosity_strength, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(mass_loc, &params.mass, SHADER_UNIFORM_FLOAT, 1);

    rlSetUniform(spiky_grad_loc, &kernel_constants.spiky_pow2_grad_scale, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(near_grad_loc, &kernel_constants.spiky_pow3_grad_scale, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(visc_loc, &kernel_constants.viscosity_scale, SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(rest_density_loc, &params.rest_density, SHADER_UNIFORM_FLOAT, 1);

    rlBindShaderBuffer(ssbo_id, 0);       // Particles (RW)
    rlBindShaderBuffer(ssbo_indices, 1);  // Indices (Read)
    rlBindShaderBuffer(ssbo_offsets, 2);  // Offsets (Read)

    force_shader.dispatch(num_groups, 1, 1);

    rlDisableShader();
}

void FluidSimulation::draw() {
    BeginShaderMode(particle_shader);
    rlEnableShader(particle_shader.id);

    Matrix mvp = MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection());
    SetShaderValueMatrix(particle_shader, particle_shader.locs[SHADER_LOC_MATRIX_MVP], mvp);

    rlBindShaderBuffer(ssbo_id, 0);

    rlEnableVertexArray(quad_mesh.vaoId);

    rlDrawVertexArrayElementsInstanced(0, quad_mesh.triangleCount * 3, 0, params.num_particles);

    rlDisableVertexArray();
    EndShaderMode();
}

void FluidSimulation::update_kernel_constants() {
    float h = params.smoothing_radius;
    float h2 = h * h;
    float h4 = h2 * h2;
    float h5 = h4 * h;
    float h8 = h4 * h4;

    kernel_constants.smoothing_radius_sq = h2;
    kernel_constants.inv_smoothing_radius = 1.0f / h;

    // Poly6: 4 / (PI * h^8)
    kernel_constants.poly6_scale = 4.0f / (PI * h8);

    // Spiky Pow3 (Near Density): 10 / (PI * h^5)
    kernel_constants.spiky_pow3_scale = 10.0f / (PI * h5);

    // Spiky Pow2 Gradient (Pressure): 12 / (PI * h^4)
    kernel_constants.spiky_pow2_grad_scale = 12.0f / (PI * h4);

    // Spiky Pow3 Gradient (Near Pressure): 30 / (PI * h^5)
    kernel_constants.spiky_pow3_grad_scale = 30.0f / (PI * h5);

    // Viscosity: 40 / (PI * h^5)
    kernel_constants.viscosity_scale = 40.0f / (PI * h5);
}

void FluidSimulation::draw_gui() {
    // GUI Panel
    int panel_width = 210;
    int panel_x = 800 - panel_width - 10;
    int panel_y = 10;
    int slider_height = 20;
    int spacing = 30;
    int y = panel_y;

    DrawRectangle(panel_x - 5, panel_y - 5, panel_width + 10, 480, Fade(DARKGRAY, 0.8f));

    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, "=== FLUID CONTROLS ===");
    y += spacing;

    // Gravity toggle
    GuiCheckBox((Rectangle){(float)panel_x, (float)y, 20, 20}, "Gravity Enabled", &params.gravity_enabled);
    y += spacing;

    // Gravity strength
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Gravity: %.1f", params.gravity));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.gravity, 0.0f, 50.0f);
    y += spacing;

    // Stiffness
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Stiffness: %.1f", params.stiffness));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.stiffness, 1.0f, 200.0f);
    y += spacing;

    // Rest Density
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Rest Density: %.1f", params.rest_density));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.rest_density, 1.0f, 100.0f);
    y += spacing;

    // Mass
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Mass: %.2f", params.mass));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.mass, 0.1f, 5.0f);
    y += spacing;

    // Viscosity Strength
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Viscosity: %.3f", params.viscosity_strength));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.viscosity_strength, 0.0f, 0.5f);
    y += spacing;

    // Smoothing Radius
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Smoothing Radius: %.2f", params.smoothing_radius));
    y += 18;
    float old_radius = params.smoothing_radius;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.smoothing_radius, 0.2f, 2.0f);
    // Update spatial hash and constants if radius changed appropriately
    if (fabs(params.smoothing_radius - old_radius) > 0.001f) {
        update_kernel_constants();
    }
    y += spacing;

    // Interaction Strength
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Interact Strength: %.1f", params.interaction_strength));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.interaction_strength, 10.0f, 5000.0f);
    y += spacing;

    // Interaction Radius
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Interact Radius: %.1f", params.interaction_radius));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.interaction_radius, 0.5f, 5.0f);
    y += spacing;

    // Near Pressure (Clavet)
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Near Pressure: %.1f", params.near_pressure_multiplier));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.near_pressure_multiplier, 0.0f, 10.0f);
    y += spacing;

    // // Boundary Damping
    // GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
    //          TextFormat("Boundary Damping: %.2f", params.boundary_damping));
    // y += 18;
    // GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
    //           &params.boundary_damping, -1.0f, 0.0f);
    // y += spacing + 10;

    // Reset button
    if (GuiButton((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, 30}, "RESET SIMULATION")) {
        reset();
    }
    y += 40;

    // Info
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Particles: %d", params.num_particles));
    y += 20;
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, "LMB: Push  RMB: Pull");
}

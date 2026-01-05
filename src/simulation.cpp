#include "simulation.h"

#include <algorithm>
#include <cmath>

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

const int WINDOW_WIDTH = 800 / PIXELS_PER_METER;
const int WINDOW_HEIGHT = 600 / PIXELS_PER_METER;

FluidSimulation::FluidSimulation() : spatial_hash(params.smoothing_radius, params.num_particles) {
    init();
    update_kernel_constants();
}

FluidSimulation::~FluidSimulation() {}

void FluidSimulation::init() { reset(); }

void FluidSimulation::reset() {
    position.clear();
    predicted_position.clear();
    velocity.clear();
    force.clear();
    density.clear();
    near_density.clear();
    pressure.clear();

    int particle_per_row = (int)sqrt(params.num_particles);
    float spacing = params.smoothing_radius * 0.9f;

    float startX = (WINDOW_WIDTH - particle_per_row * spacing) / 2.0f;
    float startY = WINDOW_HEIGHT / 4.0f;

    for (int i = 0; i < params.num_particles; i++) {
        float x = startX + (i % particle_per_row) * spacing;
        float y = startY + (i / particle_per_row) * spacing;

        x += (float)GetRandomValue(-10, 10) / 100.0f;

        x += (float)GetRandomValue(-10, 10) / 100.0f;

        position.push_back({x, y});
        predicted_position.push_back({x, y});
        velocity.push_back({0, 0});
        force.push_back({0, 0});
        density.push_back(0);
        near_density.push_back(0);
        pressure.push_back(0);
    }

    spatial_hash = SpatialHash(params.smoothing_radius, params.num_particles);
    update_kernel_constants();
}

void FluidSimulation::update(float delta_time) {
    if (delta_time > 0.03f) delta_time = 0.03f;

    int substeps = 4;
    float subtime = delta_time / (float)substeps;

    for (int i = 0; i < substeps; i++) {
        for (size_t k = 0; k < position.size(); k++) {
            // Prediction: Predict where particles will be based on current velocity
            predicted_position[k].x = position[k].x + velocity[k].x * subtime;
            predicted_position[k].y = position[k].y + velocity[k].y * subtime;
        }

        spatial_hash.update(predicted_position);

        compute_density_pressure();
        compute_forces();
        apply_mouse_force();
        integrate(subtime);
    }
}

void FluidSimulation::draw() {
    for (size_t i = 0; i < position.size(); i++) {
        Vector2 screen_pos = {position[i].x * PIXELS_PER_METER, position[i].y * PIXELS_PER_METER};

        float speed = sqrt(velocity[i].x * velocity[i].x + velocity[i].y * velocity[i].y);
        float t = speed / 8.0f;

        if (t > 1.0f) t = 1.0f;

        Color p_color = {(unsigned char)(0 + t * 150), (unsigned char)(50 + t * 205), (unsigned char)(150 + t * 105),
                         255};

        DrawCircleV(screen_pos, 4.0f, p_color);
    }
}

void FluidSimulation::apply_mouse_force() {
    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON) || IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
        Vector2 mousePos = GetMousePosition();

        if (mousePos.x > 800 - 220) return;

        float world_x = mousePos.x / PIXELS_PER_METER;
        float world_y = mousePos.y / PIXELS_PER_METER;

        // Left Click = Repel (+), Right Click = Attract (-)
        float strength =
            (IsMouseButtonDown(MOUSE_LEFT_BUTTON) ? params.interaction_strength : -params.interaction_strength);

        for (size_t i = 0; i < position.size(); i++) {
            float dx = world_x - position[i].x;
            float dy = world_y - position[i].y;
            float distSq = dx * dx + dy * dy;

            if (distSq < params.interaction_radius * params.interaction_radius) {
                float dist = sqrt(distSq);
                if (dist < 0.001f) dist = 0.001f;

                float dirX = dx / dist;
                float dirY = dy / dist;

                float factor = 1.0f - (dist / params.interaction_radius);

                // F = Strength * Factor
                // Acceleration = F / Density (approx mass/density)
                float force_val = strength * factor;

                force[i].x -= dirX * force_val * density[i];
                force[i].y -= dirY * force_val * density[i];
            }
        }
    }
}

void FluidSimulation::compute_density_pressure() {
    for (size_t i = 0; i < position.size(); i++) {
        density[i] = 0.0f;
        near_density[i] = 0.0f;

        spatial_hash.for_each_neighbor(predicted_position[i], [&](int j) {
            float dx = predicted_position[j].x - predicted_position[i].x;
            float dy = predicted_position[j].y - predicted_position[i].y;

            float distance_squared = dx * dx + dy * dy;

            if (distance_squared < kernel_constants.smoothing_radius_sq) {
                float distance = sqrt(distance_squared);
                density[i] += params.mass * poly6_kernel(distance_squared);
                near_density[i] += params.mass * spiky_pow3_kernel(distance);
            }
        });

        if (density[i] < 0.0001f) density[i] = 0.0001f;

        pressure[i] = params.stiffness * (density[i] - params.rest_density);

        if (pressure[i] < 0.0f) pressure[i] = 0.0f;
    }
}

void FluidSimulation::compute_forces() {
    for (size_t i = 0; i < position.size(); i++) {
        Vector2 pressure_force = {0.0f, 0.0f};
        Vector2 viscosity_force = {0.0f, 0.0f};

        spatial_hash.for_each_neighbor(predicted_position[i], [&](int j) {
            if (i == (size_t)j) return;

            float dx = predicted_position[j].x - predicted_position[i].x;
            float dy = predicted_position[j].y - predicted_position[i].y;
            float distance_squared = dx * dx + dy * dy;

            if (distance_squared > 0 && distance_squared < kernel_constants.smoothing_radius_sq) {
                float distance = sqrt(distance_squared);
                float direction_x = dx / distance;
                float direction_y = dy / distance;

                // --- Pressure Forces (Double Density Relaxation) ---
                float shared_pressure = (pressure[i] + pressure[j]) / 2.0f;
                float shared_near_pressure = (near_density[i] * params.near_pressure_multiplier +
                                              near_density[j] * params.near_pressure_multiplier) /
                                             2.0f;

                // Regular Pressure: Uses Linear Gradient (spiky_pow2_gradient)
                float slope = spiky_pow2_gradient(distance);  // Linear
                // Near Pressure: Uses Quadratic Gradient (spiky_pow3_gradient)
                float near_slope = spiky_pow3_gradient(distance);  // Quadratic

                // F = -m * [ (P+P)/2 * Grad + (P_near+P_near)/2 * Grad_Near ] / rho
                float p_term = shared_pressure * slope;
                float near_p_term = shared_near_pressure * near_slope;

                float term = params.mass * (p_term + near_p_term) / density[j];

                pressure_force.x += term * direction_x;
                pressure_force.y += term * direction_y;

                // --- Viscosity Force ---
                // Fv = mu * m * (vj - vi) / rho_j * Laplacian(W)
                float visc_laplacian = viscosity_kernel(distance);
                float vx_diff = velocity[j].x - velocity[i].x;
                float vy_diff = velocity[j].y - velocity[i].y;

                float visc_term = params.viscosity_strength * params.mass * visc_laplacian / density[j];
                viscosity_force.x += vx_diff * visc_term;
                viscosity_force.y += vy_diff * visc_term;
            }
        });

        force[i].x = pressure_force.x + viscosity_force.x;
        force[i].y = pressure_force.y + viscosity_force.y;
    }
}

void FluidSimulation::integrate(float delta_time) {
    for (size_t i = 0; i < position.size(); i++) {
        // Apply gravity if enabled
        if (params.gravity_enabled) {
            force[i].y += params.gravity * density[i];
        }

        Vector2 acceleration = {force[i].x / density[i], force[i].y / density[i]};

        velocity[i].x += acceleration.x * delta_time;
        velocity[i].y += acceleration.y * delta_time;

        position[i].x += velocity[i].x * delta_time;
        position[i].y += velocity[i].y * delta_time;

        float boundary_margin = 0.1f;

        if (position[i].x > WINDOW_WIDTH - boundary_margin) {
            position[i].x = WINDOW_WIDTH - boundary_margin;
            velocity[i].x *= params.boundary_damping;
        }
        if (position[i].x < boundary_margin) {
            position[i].x = boundary_margin;
            velocity[i].x *= params.boundary_damping;
        }
        if (position[i].y > WINDOW_HEIGHT - boundary_margin) {
            position[i].y = WINDOW_HEIGHT - boundary_margin;
            velocity[i].y *= params.boundary_damping;
        }
        if (position[i].y < boundary_margin) {
            position[i].y = boundary_margin;
            velocity[i].y *= params.boundary_damping;
        }
    }
}

float FluidSimulation::poly6_kernel(float distance_squared) {
    if (distance_squared >= kernel_constants.smoothing_radius_sq) {
        return 0.0f;
    }
    float diff = kernel_constants.smoothing_radius_sq - distance_squared;

    return kernel_constants.poly6_scale * diff * diff * diff;
}

float FluidSimulation::spiky_pow2_gradient(float distance) {
    if (distance >= params.smoothing_radius || distance <= 0.0f) {
        return 0.0f;
    }
    float diff = params.smoothing_radius - distance;

    return -12.0f / (PI * pow(params.smoothing_radius, 4)) * diff;
}

float FluidSimulation::spiky_pow3_kernel(float distance) {
    if (distance >= params.smoothing_radius) return 0.0f;
    float diff = params.smoothing_radius - distance;

    return kernel_constants.spiky_pow3_scale * diff * diff * diff;
}

float FluidSimulation::spiky_pow3_gradient(float distance) {
    if (distance >= params.smoothing_radius || distance <= 0.0f) return 0.0f;
    float diff = params.smoothing_radius - distance;

    return kernel_constants.spiky_pow3_grad_scale * diff * diff;
}

float FluidSimulation::viscosity_kernel(float distance) {
    if (distance >= params.smoothing_radius || distance < 0.0f) {
        return 0.0f;
    }
    float diff = params.smoothing_radius - distance;

    return kernel_constants.viscosity_scale * diff;
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

    // Spiky Pow2 Gradient (Pressure): -12 / (PI * h^4)
    kernel_constants.spiky_pow2_grad_scale = -12.0f / (PI * h4);

    // Spiky Pow3 Gradient (Near Pressure): -30 / (PI * h^5)
    kernel_constants.spiky_pow3_grad_scale = -30.0f / (PI * h5);

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
        if (fabs(params.smoothing_radius - old_radius) > 0.01f) {
            spatial_hash = SpatialHash(params.smoothing_radius, params.num_particles);
        }
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

    // Boundary Damping
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Boundary Damping: %.2f", params.boundary_damping));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "",
              &params.boundary_damping, -1.0f, 0.0f);
    y += spacing + 10;

    // Reset button
    if (GuiButton((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, 30}, "RESET SIMULATION")) {
        reset();
    }
    y += 40;

    // Info
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20},
             TextFormat("Particles: %d", (int)position.size()));
    y += 20;
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, "LMB: Push  RMB: Pull");
}

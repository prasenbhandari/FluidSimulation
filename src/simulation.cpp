#include "simulation.h"
#include <cmath>
#include <algorithm>

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

const int WINDOW_WIDTH = 800 / PIXELS_PER_METER;
const int WINDOW_HEIGHT = 600 / PIXELS_PER_METER;

const float PREDICTION_FACTOR = 1.0f / 120.0f;

FluidSimulation::FluidSimulation(): spatial_hash(params.smoothing_radius, params.num_particles) { 
    init(); 
}


FluidSimulation::~FluidSimulation() {}

void FluidSimulation::init() { reset(); }

void FluidSimulation::reset() {
    particles.clear();

    int particle_per_row = (int)sqrt(params.num_particles);
    float spacing = params.smoothing_radius * 0.9f;

    float startX = (WINDOW_WIDTH - particle_per_row * spacing) / 2.0f;
    float startY = WINDOW_HEIGHT / 4.0f;

    for (int i = 0; i < params.num_particles; i++){
        float x = startX + (i % particle_per_row) * spacing;
        float y = startY + (i / particle_per_row) * spacing;

        x += (float)GetRandomValue(-10, 10) / 100.0f;

        particles.emplace_back(x, y);
    }
    
    // Recreate spatial hash with new parameters
    spatial_hash = SpatialHash(params.smoothing_radius, params.num_particles);
}


void FluidSimulation::update(float delta_time){
    handle_input();

    if (delta_time > 0.03f) delta_time = 0.03f;

    int substeps = 4;
    float subtime = delta_time / (float)substeps;

    for (int i = 0; i < substeps; i++) {

        for (auto& p: particles){
            p.predicted_position.x = p.position.x + p.velocity.x * PREDICTION_FACTOR;
            p.predicted_position.y = p.position.y + p.velocity.y * PREDICTION_FACTOR;
        }

        spatial_hash.update(particles);

        compute_density_pressure();
        compute_forces();
        integrate(subtime);
    }
}


void FluidSimulation::draw() {
    for (const auto& p: particles){
        Vector2 screen_pos = {
            p.position.x * PIXELS_PER_METER,
            p.position.y * PIXELS_PER_METER
        };

        float speed = sqrt(p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y);
        float t = speed / 8.0f;

        if (t > 1.0f) t = 1.0f;

        Color p_color = {
            (unsigned char)(0 + t * 150),
            (unsigned char)(50 + t * 205),
            (unsigned char)(150 + t * 105),
            255
        };

        DrawCircleV(screen_pos, 4.0f, p_color);        
    }
}


void FluidSimulation::handle_input() {
    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON) || IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
        Vector2 mousePos = GetMousePosition();
        
        // Don't interact if mouse is over GUI panel (right side)
        if (mousePos.x > 800 - 220) return;

        float world_x = mousePos.x / PIXELS_PER_METER;
        float world_y = mousePos.y / PIXELS_PER_METER;

        // Left Click = Repel (+), Right Click = Attract (-)
        float strength = IsMouseButtonDown(MOUSE_LEFT_BUTTON) ? params.interaction_strength : -params.interaction_strength;

        for (auto& p : particles) {
            float dx = world_x - p.position.x;
            float dy = world_y - p.position.y;
            float distSq = dx*dx + dy*dy;

            if (distSq < params.interaction_radius * params.interaction_radius) {
                float dist = sqrt(distSq);
                if (dist < 0.001f) dist = 0.001f; 

                float dirX = dx / dist;
                float dirY = dy / dist;

                float factor = 1.0f - (dist / params.interaction_radius);

                p.velocity.x -= dirX * strength * factor;
                p.velocity.y -= dirY * strength * factor;
            }
        }
    }
}


void FluidSimulation::compute_density_pressure() {
    for (size_t i = 0; i < particles.size(); i++){
        Particle& pi = particles[i];
        pi.density = 0.0f;

        spatial_hash.for_each_neighbor(pi.predicted_position, [&](int j){
            const Particle& pj = particles[j];

            float dx = pj.predicted_position.x - pi.predicted_position.x;
            float dy = pj.predicted_position.y - pi.predicted_position.y;

            float distance_squared = dx * dx + dy * dy;

            if(distance_squared < params.smoothing_radius * params.smoothing_radius){
                pi.density += params.mass * poly6_kernel(distance_squared);
            }
        });

        if (pi.density < 0.0001f) pi.density = 0.0001f;
        // Use positive pressure only - particles repel when density > 0
        // This creates gas-like spreading behavior
        pi.pressure = params.stiffness * pi.density;
    }
}


void FluidSimulation::compute_forces() {
    for(size_t i = 0; i < particles.size(); i++){
        Particle& pi = particles[i];
        Vector2 pressure_force = {0.0f, 0.0f};

        spatial_hash.for_each_neighbor(pi.predicted_position, [&](int j){
            if (i == (size_t)j) return;

            const Particle& pj = particles[j];

            float dx = pj.predicted_position.x - pi.predicted_position.x;
            float dy = pj.predicted_position.y - pi.predicted_position.y;
            float distance_squared = dx * dx + dy * dy;

            if (distance_squared > 0 && distance_squared < params.smoothing_radius * params.smoothing_radius){
                float distance = sqrt(distance_squared);
                float direction_x = dx / distance;
                float direction_y = dy / distance;

                // Average pressure, divide by both densities
                float shared_pressure = (pi.pressure + pj.pressure) / 2.0f;
                float slope = spiky_gradient(distance);

                pressure_force.x += params.mass * shared_pressure * slope * direction_x / pj.density;
                pressure_force.y += params.mass * shared_pressure * slope * direction_y / pj.density;
            }

        });

        pi.force = pressure_force;
    }
}


void FluidSimulation::integrate(float delta_time) {
    for (auto& p : particles) {
        // Apply gravity if enabled
        if (params.gravity_enabled) {
            p.force.y += params.gravity * p.density;
        }

        Vector2 acceleration = {
            p.force.x / p.density,
            p.force.y / p.density
        };

        p.velocity.x += acceleration.x * delta_time;
        p.velocity.y += acceleration.y * delta_time;

        p.position.x += p.velocity.x * delta_time;
        p.position.y += p.velocity.y * delta_time;

        float boundary_margin = 0.1f;

        if (p.position.x > WINDOW_WIDTH - boundary_margin) {
            p.position.x = WINDOW_WIDTH - boundary_margin;
            p.velocity.x *= params.boundary_damping;
        }
        if (p.position.x < boundary_margin) {
            p.position.x = boundary_margin;
            p.velocity.x *= params.boundary_damping;
        }
        if (p.position.y > WINDOW_HEIGHT - boundary_margin) {
            p.position.y = WINDOW_HEIGHT - boundary_margin;
            p.velocity.y *= params.boundary_damping;
        }
        if (p.position.y < boundary_margin) {
            p.position.y = boundary_margin;
            p.velocity.y *= params.boundary_damping;
        }
    }
}


float FluidSimulation::poly6_kernel(float distance_squared) {
    if (distance_squared >= params.smoothing_radius * params.smoothing_radius) {
        return 0.0f;
    }
    float h2 = params.smoothing_radius * params.smoothing_radius;
    float diff = h2 - distance_squared;

    // 2D Poly6 kernel: 4 / (PI * h^8) * (h^2 - r^2)^3
    return (4.0f / (PI * pow(params.smoothing_radius, 8))) * pow(diff, 3);
}

float FluidSimulation::spiky_gradient(float distance) {
    if (distance >= params.smoothing_radius || distance == 0.0f) {
        return 0.0f;
    }
    float diff = params.smoothing_radius - distance;

    // 2D Spiky derivative: -12 / (PI * h^4) * (h - r)
    return -12.0f / (PI * pow(params.smoothing_radius, 4)) * diff;
}

void FluidSimulation::draw_gui() {
    // GUI Panel on the right side
    int panel_width = 210;
    int panel_x = 800 - panel_width - 10;
    int panel_y = 10;
    int slider_height = 20;
    int spacing = 30;
    int y = panel_y;
    
    // Semi-transparent background
    DrawRectangle(panel_x - 5, panel_y - 5, panel_width + 10, 480, Fade(DARKGRAY, 0.8f));
    
    // Title
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, "=== FLUID CONTROLS ===");
    y += spacing;
    
    // Gravity toggle
    GuiCheckBox((Rectangle){(float)panel_x, (float)y, 20, 20}, "Gravity Enabled", &params.gravity_enabled);
    y += spacing;
    
    // Gravity strength
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Gravity: %.1f", params.gravity));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "", &params.gravity, 0.0f, 50.0f);
    y += spacing;
    
    // Stiffness
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Stiffness: %.1f", params.stiffness));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "", &params.stiffness, 1.0f, 200.0f);
    y += spacing;
    
    // Rest Density
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Rest Density: %.1f", params.rest_density));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "", &params.rest_density, 1.0f, 100.0f);
    y += spacing;
    
    // Mass
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Mass: %.2f", params.mass));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "", &params.mass, 0.1f, 5.0f);
    y += spacing;
    
    // Smoothing Radius
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Smoothing Radius: %.2f", params.smoothing_radius));
    y += 18;
    float old_radius = params.smoothing_radius;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "", &params.smoothing_radius, 0.2f, 2.0f);
    // Update spatial hash if radius changed significantly
    if (fabs(params.smoothing_radius - old_radius) > 0.01f) {
        spatial_hash = SpatialHash(params.smoothing_radius, params.num_particles);
    }
    y += spacing;
    
    // Interaction Strength
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Interact Strength: %.1f", params.interaction_strength));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "", &params.interaction_strength, 10.0f, 200.0f);
    y += spacing;
    
    // Interaction Radius
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Interact Radius: %.1f", params.interaction_radius));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "", &params.interaction_radius, 0.5f, 5.0f);
    y += spacing;
    
    // Boundary Damping
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Boundary Damping: %.2f", params.boundary_damping));
    y += 18;
    GuiSlider((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, (float)slider_height}, "", "", &params.boundary_damping, -1.0f, 0.0f);
    y += spacing + 10;
    
    // Reset button
    if (GuiButton((Rectangle){(float)panel_x, (float)y, (float)panel_width - 10, 30}, "RESET SIMULATION")) {
        reset();
    }
    y += 40;
    
    // Info
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, TextFormat("Particles: %d", (int)particles.size()));
    y += 20;
    GuiLabel((Rectangle){(float)panel_x, (float)y, (float)panel_width, 20}, "LMB: Push  RMB: Pull");
}

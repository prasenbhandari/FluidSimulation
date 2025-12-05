#include "simulation.h"
#include <cmath>
#include <algorithm>


const int NUM_PARTICLES = 800;
const float SMOOTHING_RADIUS = 35.0f;
const float GRAVITY = 9.84f;
const float MASS = 1.0f;
const float REST_DENSITY = 0.005f;
const float DAMPING = 0.97f;
const float BOUNDRY_DAMPING = -0.5f;
const float STIFFNESS = 1200.0f;

const int INTERACTION_RADIUS = 100;
const float INTERACTION_STRENGTH = 20.0f;

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;

const float PREDICTION_FACTOR = 1.0f / 120.0f;

FluidSimulation::FluidSimulation() { init(); }
FluidSimulation::~FluidSimulation() {}

void FluidSimulation::init() { reset(); }

void FluidSimulation::reset() {
    particles.clear();

    int particle_per_row = (int)sqrt(NUM_PARTICLES);
    float spacing = SMOOTHING_RADIUS * 0.5f;

    float startX = (WINDOW_WIDTH - particle_per_row * spacing) / 2.0f;
    float startY = WINDOW_HEIGHT / 4.0f;

    for (int i = 0; i < NUM_PARTICLES; i++){
        float x = startX + (i % particle_per_row) * spacing;
        float y = startY + (i / particle_per_row) * spacing;

        x += (float)GetRandomValue(-10, 10) / 10.f;

        particles.emplace_back(x, y);
    }
}


void FluidSimulation::update(float delta_time){
    int substeps = 2;
    float subtime = delta_time / (float)substeps;

    for (int i = 0; i < substeps; i++) {

        for (auto& p: particles){
            p.predicted_position.x = p.position.x + p.velocity.x * subtime;
            p.predicted_position.y = p.position.y + p.velocity.y * subtime;
        }

        compute_density_pressure();
        compute_forces();
        integrate(subtime);
    }
}


void FluidSimulation::draw() {
    for (const auto& p: particles){
        float speed = sqrt(p.velocity.x * p.velocity.x + p.velocity.y * p.velocity.y);
        float t = speed / 200.0f;

        if (t > 1.0f) t = 1.0f;

        Color p_color = {
            (unsigned char)(0 + t * 150),
            (unsigned char)(50 + t * 205),
            (unsigned char)(150 + t * 105),
            255
        };

        DrawCircleV(p.position, 4.0f, p_color);        
    }
}


void FluidSimulation::handle_input() {
    if (IsMouseButtonDown(MOUSE_LEFT_BUTTON) || IsMouseButtonDown(MOUSE_RIGHT_BUTTON)) {
        Vector2 mousePos = GetMousePosition();
        // Left Click = Repel (+), Right Click = Attract (-)
        float strength = IsMouseButtonDown(MOUSE_LEFT_BUTTON) ? INTERACTION_STRENGTH : -INTERACTION_STRENGTH;

        for (auto& p : particles) {
            float dx = mousePos.x - p.position.x;
            float dy = mousePos.y - p.position.y;
            float distSq = dx*dx + dy*dy;

            if (distSq < INTERACTION_RADIUS * INTERACTION_RADIUS) {
                float dist = sqrt(distSq);
                if (dist < 0.01f) dist = 0.01f; 

                float dirX = dx / dist;
                float dirY = dy / dist;

                float factor = 1.0f - (dist / INTERACTION_RADIUS);

                p.velocity.x -= dirX * strength * factor;
                p.velocity.y -= dirY * strength * factor;
            }
        }
    }
}


void FluidSimulation::compute_density_pressure() {
    for (auto& pi : particles) {
        pi.density = 0.0f;

        for (const auto& pj : particles) {
            float dx = pj.predicted_position.x - pi.predicted_position.x;
            float dy = pj.predicted_position.y - pi.predicted_position.y;
            float distance_squared = dx * dx + dy * dy;

            if (distance_squared < SMOOTHING_RADIUS * SMOOTHING_RADIUS) {
                pi.density += MASS * poly6_kernel(distance_squared);
            }
        }

        if (pi.density < 0.0001f) pi.density = 0.0001f;

        pi.pressure = STIFFNESS * (pi.density - REST_DENSITY);
    }
}


void FluidSimulation::compute_forces() {
    for(size_t i = 0; i < particles.size(); i++){
        Particle& pi = particles[i];
        Vector2 pressure_force = {0.0f, 0.0f};

        for (size_t j = 0; j < particles.size(); j++){
            if (i == j) continue;

            const Particle& pj = particles[j];
            float dx = pj.predicted_position.x - pi.predicted_position.x;
            float dy = pj.predicted_position.y - pi.predicted_position.y;

            float distance_squared = dx * dx + dy * dy;

            if (distance_squared > 0 && distance_squared < SMOOTHING_RADIUS * SMOOTHING_RADIUS) {
                float distance = sqrt(distance_squared);
                float direction_x = dx / sqrt(distance_squared);
                float direction_y = dy / sqrt(distance_squared);

                float shared_pressure = (pi.pressure + pj.pressure) / (2.0f * pj.density);
                float slope = spiky_gradient(distance);

                pressure_force.x += -MASS * shared_pressure * slope * direction_x;
                pressure_force.y += -MASS * shared_pressure * slope * direction_y;
            }
        }

        pi.force = pressure_force;
    }
}


void FluidSimulation::integrate(float delta_time) {
    for (auto& p : particles) {
        // Apply gravity
        p.force.y += GRAVITY * p.density;

        
        Vector2 acceleration = {
            p.force.x / p.density,
            p.force.y / p.density
        };

        p.velocity.x += acceleration.x * delta_time;
        p.velocity.y += acceleration.y * delta_time;

        p.position.x += p.velocity.x * delta_time;
        p.position.y += p.velocity.y * delta_time;

        if (p.position.x > WINDOW_WIDTH - 10.0f) {
            p.position.x = WINDOW_WIDTH - 10.0f;
            p.velocity.x *= BOUNDRY_DAMPING;
        }
        if (p.position.x < 10.0f) {
            p.position.x = 10.0f;
            p.velocity.x *= BOUNDRY_DAMPING;
        }
        if (p.position.y > WINDOW_HEIGHT - 10.0f) {
            p.position.y = WINDOW_HEIGHT - 10.0f;
            p.velocity.y *= BOUNDRY_DAMPING;
        }
        if (p.position.y < 10.0f) {
            p.position.y = 10.0f;
            p.velocity.y *= BOUNDRY_DAMPING;
        }
    }
}


float FluidSimulation::poly6_kernel(float distance_squared) {
    if (distance_squared >= SMOOTHING_RADIUS * SMOOTHING_RADIUS) {
        return 0.0f;
    }
    float h2 = SMOOTHING_RADIUS * SMOOTHING_RADIUS;
    float diff = h2 - distance_squared;

    return (315.0f / (64.0f * PI * pow(SMOOTHING_RADIUS, 9))) * pow(diff, 3);
}

float FluidSimulation::spiky_gradient(float distance) {
    if (distance >= SMOOTHING_RADIUS || distance == 0.0f) {
        return 0.0f;
    }
    float diff = SMOOTHING_RADIUS - distance;

    return -45.0f / (PI * pow(SMOOTHING_RADIUS, 6)) * diff * diff;
}

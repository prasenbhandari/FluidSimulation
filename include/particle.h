#pragma once

#include "raylib.h"


struct Particle {
    Vector2 position;
    Vector2 predicted_position;
    Vector2 velocity;
    Vector2 force;
    float density;
    float pressure;

    Particle(float x, float y) {
        position = {x, y};
        predicted_position = {x, y};
        velocity = {0.0f, 0.0f};
        density = 0.0f;
        pressure = 0.0f;
    }
};

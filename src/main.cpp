
#include "raylib.h"
#include "rlgl.h"
#include "simulation.h"

const int WINDOW_WIDTH = WINDOW_WIDTH_PIXELS;
const int WINDOW_HEIGHT = WINDOW_HEIGHT_PIXELS;

int main() {
    SetConfigFlags(FLAG_VSYNC_HINT);
    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Fluid Simulation");
    SetTargetFPS(60);

    FluidSimulation sim;

    Camera2D camera = {0};
    camera.zoom = PIXELS_PER_METER;

    const float FIXED_DT = 1.0f / 60.0f;

    while (!WindowShouldClose()) {
        sim.update(FIXED_DT);

        BeginDrawing();
        ClearBackground(BLACK);

        BeginMode2D(camera);
        sim.draw();
        EndMode2D();

        DrawFPS(10, 10);
        sim.draw_gui();

        EndDrawing();
    }

    CloseWindow();

    return 0;
}

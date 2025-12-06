#include "simulation.h"
#include "raylib.h"

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;


int main(){
    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Fluid Simulation");

    FluidSimulation sim;
    sim.init();

    SetTargetFPS(60);

    while (!WindowShouldClose()){
        float delta_time = GetFrameTime();

        if (delta_time > 0.05f) delta_time = 0.05f;

        sim.update(delta_time);

        BeginDrawing();
        ClearBackground(BLACK);
        DrawFPS(10, 10);
        sim.draw();
        sim.draw_gui();  // Draw raygui controls

        EndDrawing();
    }

    CloseWindow();

    return 0;
}

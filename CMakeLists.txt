cmake_minimum_required(VERSION 3.14)

project(FluidSimRaylib VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    raylib
    GIT_REPOSITORY https://github.com/raysan5/raylib.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)

# Force OpenGL 4.3 for Raylib
set(OPENGL_VERSION "4.3" CACHE STRING "OpenGL 4.3" FORCE)

FetchContent_MakeAvailable(raylib)

include_directories(include)

set(SOURCES
    src/main.cpp
    src/compute_shader.cpp
    src/simulation.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PRIVATE raylib)
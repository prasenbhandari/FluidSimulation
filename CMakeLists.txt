cmake_minimum_required(VERSION 3.10)

project(FluidSimRaylib VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Dependencies ---
# Attempt to find the Raylib package installed on the system
find_package(raylib REQUIRED)

include_directories(include)

set(SOURCES
    src/main.cpp
    src/simulation.cpp
    src/particle.cpp
)

# --- Build Target ---
add_executable(${PROJECT_NAME} ${SOURCES})

# --- Linking ---
target_link_libraries(${PROJECT_NAME} PRIVATE raylib)

# Specific linking for Windows/Linux/macOS if needed
if(WIN32)
    target_link_libraries(${PROJECT_NAME} PRIVATE gdi32 winmm)
elseif(UNIX AND NOT APPLE)
    target_link_libraries(${PROJECT_NAME} PRIVATE m pthread dl rt X11)
endif()


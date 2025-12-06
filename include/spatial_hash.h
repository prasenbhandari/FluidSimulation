#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include "particle.h"

struct spatial_entry {
    int particle_index;
    unsigned int cell_key;
};

class SpatialHash {
public:
    SpatialHash(float spacing, int table_size);

    void update(const std::vector<Particle>& particles);

    template<typename Func>
    void for_each_neighbor(Vector2 position, Func func){
        int cell_x = (int)std::floor(position.x / spacing);
        int cell_y = (int)std::floor(position.y / spacing);

        for (int x = -1; x <= 1; x++){
            for(int y = -1; y <= 1; y++){
                unsigned int key = hash(cell_x + x, cell_y + y);
                int start = start_indices[key];

                if (start == -1) continue;

                for (size_t i = start; i < num_particles && spatial_lookup[i].cell_key == key; i++){
                    func(spatial_lookup[i].particle_index);
                }
            }
        }
    }

private:
    float spacing;
    int table_size;
    size_t num_particles = 0;

    std::vector<spatial_entry> spatial_lookup;

    std::vector<int> start_indices;

    unsigned int hash(int x, int y);
};
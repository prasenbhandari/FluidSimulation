#include "spatial_hash.h"

const unsigned int PRIME_1 = 15823;
const unsigned int PRIME_2 = 9737333;

SpatialHash::SpatialHash(float spacing, int table_size): 
    spacing(spacing), table_size(table_size) {
    start_indices.resize(table_size, -1);
}


unsigned int SpatialHash::hash(int x, int y) {
    unsigned int hash_value = (unsigned int)(x * PRIME_1) ^ (unsigned int)(y * PRIME_2);
    return hash_value % table_size;
}


void SpatialHash::update(const std::vector<Particle>& particles){
    std::fill(start_indices.begin(), start_indices.end(), -1);

    num_particles = particles.size();
    
    if(spatial_lookup.size() < num_particles){
        spatial_lookup.resize(num_particles);
    }

    for (size_t i = 0; i < num_particles; i++){
        int cell_x = (int)std::floor(particles[i].predicted_position.x / spacing);
        int cell_y = (int)std::floor(particles[i].predicted_position.y / spacing);

        spatial_lookup[i].particle_index = (int)i;
        spatial_lookup[i].cell_key = hash(cell_x, cell_y);
    }

    std::sort(spatial_lookup.begin(), spatial_lookup.begin() + num_particles,
        [](const spatial_entry& a, const spatial_entry& b){
            return a.cell_key < b.cell_key;
        }
    );


    for(size_t i = 0; i < num_particles; i++){
        unsigned int key = spatial_lookup[i].cell_key;

        unsigned int prev_key = (i == 0) ? (unsigned int)(-1) : spatial_lookup[i - 1].cell_key;

        if (key != prev_key){
            start_indices[key] = (int)i;
        }
    }
}
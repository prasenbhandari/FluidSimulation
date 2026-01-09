#pragma once

#include <string>

#include "raylib.h"

class ComputeShader {
public:
    ComputeShader();
    ~ComputeShader();

    bool load(const std::string& start_file_path);

    void dispatch(unsigned int x, unsigned int y, unsigned int z);

    static void memory_barrier(unsigned int barriers);

    unsigned int get_program_id() const;
    unsigned int get_shader_id() const;

private:
    unsigned int program_id;
    unsigned int shader_id;
};

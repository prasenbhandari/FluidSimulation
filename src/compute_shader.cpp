#include "compute_shader.h"

#include <raylib.h>

#include "external/glad.h"
#include "rlgl.h"

ComputeShader::ComputeShader() : program_id(0), shader_id(0) {}

ComputeShader::~ComputeShader() { rlUnloadShaderProgram(program_id); }

bool ComputeShader::load(const std::string& file_path) {
    char* code = LoadFileText(file_path.c_str());
    if (code == nullptr) {
        return false;
    }

    shader_id = rlCompileShader(code, RL_COMPUTE_SHADER);

    program_id = rlLoadComputeShaderProgram(shader_id);

    UnloadFileText(code);

    return program_id != 0;
}

void ComputeShader::dispatch(unsigned int x, unsigned int y, unsigned int z) {
    rlComputeShaderDispatch(x, y, z);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}

unsigned int ComputeShader::get_program_id() const { return program_id; }

unsigned int ComputeShader::get_shader_id() const { return shader_id; }

#ifndef SHADER_HPP
#define SHADER_HPP

#include <string>

#include <GLFW/glfw3.h>

GLuint loadShader(std::string vertex_fragment_file_path, std::string macro_prepends, bool geometry);
void loadBasicShader(std::string path);
void deleteShaderPrograms();

namespace shader {
    enum class TYPE {
        BASIC_SHADER = 0,
        NUM_SHADER_TYPES,
    };

    extern GLuint basic_program;
    extern struct BasicUniforms {
        GLuint mvp, model, sun_color, sun_direction, camera_position, albedo;
    } basic_uniforms;

    extern std::string glsl_version;
}

#endif

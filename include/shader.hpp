#ifndef SHADER_HPP
#define SHADER_HPP

#include <string>

#include <GLFW/glfw3.h>

GLuint loadShader(std::string_view vertex_fragment_file_path, std::string_view macro_prepends, bool geometry);
void loadBasicShader(std::string_view path);
void loadBasicInstancedShader(std::string_view path);
void loadNullShader(std::string_view path);
void loadNullInstancedShader(std::string_view path);
void deleteShaderPrograms();

namespace shader {
    enum class TYPE {
        BASIC_SHADER = 0,
        BASIC_INSTANCED_SHADER = 1,
        NUM_SHADER_TYPES,
    };

    extern GLuint basic_program;
    extern struct BasicUniforms {
        GLuint mvp, model, sun_color, sun_direction, camera_position, albedo;
    } basic_uniforms;

    extern GLuint basic_instanced_program;
    extern struct BasicUniforms basic_instanced_uniforms;

    extern GLuint null_program;
    extern struct NullUniforms {
        GLuint mvp, line_width;
    } null_uniforms;

    extern GLuint null_instanced_program;
    extern struct NullUniforms null_instanced_uniforms;

    extern std::string glsl_version;
}

#endif

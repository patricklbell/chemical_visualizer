#ifndef SHADER_CORE_HPP
#define SHADER_CORE_HPP

#include <unordered_map>
#include <filesystem>

#ifdef EMSCRIPTEN
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#else
#include <GL/glew.h>
#endif

#include <glm/glm.hpp>

struct Shader {
    ~Shader();
    bool load_file(std::string_view path);
    bool load_file_compile(std::string_view path, std::string_view _prepend = "");
    bool ready();
    bool compile(std::string_view _prepend = "");
    void clear();
    bool live_update();
    bool apply_macros();    // Returns true if recompile was necessary
    bool set_macro(std::string_view macro, bool value, bool activate = true);

    constexpr GLuint program() { return active_program; }
    void bind();

    GLuint uniform(std::string_view name);
    void uniform(std::string_view name, const glm::mat4 &v);
    void uniform(std::string_view name, const glm::mat3 &v);
    void uniform(std::string_view name, const glm::vec4 &v);
    void uniform(std::string_view name, const glm::vec3 &v);
    void uniform(std::string_view name, const glm::vec2 &v);

    std::string handle = "";

    enum class Type : uint64_t {
        VERTEX   = GL_VERTEX_SHADER,
        FRAGMENT = GL_FRAGMENT_SHADER,
        // Emscripten doesn't seem to support these shader types
        // GEOMETRY                = GL_GEOMETRY_SHADER,
        // COMPUTE                 = GL_COMPUTE_SHADER,
        // TESSELLATION_CONTROL    = GL_TESS_CONTROL_SHADER,
        // TESSELLATION_EVALUATION = GL_TESS_EVALUATION_SHADER,
        NONE = GL_FALSE,
    };
    struct FileDependency {
        std::string path;
        std::filesystem::file_time_type last_write_time;
    };

   private:
    bool file_loaded = false;
    std::unordered_map<std::string, FileDependency> dependencies;
    std::unordered_map<Type, std::string> type_to_chunk;

    std::unordered_map<uint64_t, std::unordered_map<std::string, GLuint>> programs_uniforms;
    std::string prepend = "";
    std::unordered_map<std::string, bool> macros;
    std::unordered_map<uint64_t, GLuint> programs;

    uint64_t active_hash  = 0;
    GLuint active_program = GL_FALSE;
};

#endif    // SHADER_CORE_HPP
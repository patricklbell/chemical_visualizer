#include <glm/glm.hpp>
#include <iostream>
#include <shader/globals.hpp>

#include "graphics.hpp"

namespace Shaders {
    Shader basic;
    Shader basic_instanced;
};    // namespace Shaders

using namespace Shaders;

[[nodiscard]] bool init_shaders() {
#ifdef EMSCRIPTEN
    std::string prepend =
        "#version 300 es\n"
        "#define GLES\n"
        "#ifdef GL_FRAGMENT_PRECISION_HIGH\n"
        "precision highp float;\n"
        "#else\n"
        "precision mediump float;\n"
        "#endif\n\0";
#else
    std::string prepend = "#version 330 core\n\0";
#endif

    bool success = basic.load_file_compile("data/shaders/basic.gl", prepend);
    success &= basic_instanced.load_file_compile("data/shaders/basic_instanced.gl", prepend);

    return success;
}

[[nodiscard]] bool live_update_shaders() {
    bool success = true;

    success &= basic.live_update();
    success &= basic_instanced.live_update();

    return success;
}
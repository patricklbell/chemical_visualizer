#ifndef SHADER_GLOBALS_HPP
#define SHADER_GLOBALS_HPP

#include <shader/core.hpp>

namespace Shaders {
    extern Shader basic;
    extern Shader basic_instanced;
};    // namespace Shaders

bool init_shaders();
bool live_update_shaders();

#endif    // SHADER_GLOBALS_HPP
#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <vector>

#ifdef EMSCRIPTEN
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#else
#include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>

#include <camera/core.hpp>
#include <entities.hpp>
#include <graphics.hpp>
#include <loader.hpp>

namespace Globals {
    extern GLFWwindow* window;

    extern Camera camera;
    extern Entities entities;
    extern DrawSettings draw_settings;

    extern MolFile mol;

    extern PdbDictionary pdb_dictionary;
    extern PdbDrawSettings pdb_settings;
    extern int pdb_selected_model;
    extern PdbFile pdb;
}    // namespace Globals

#endif    // GLOBALS_HPP

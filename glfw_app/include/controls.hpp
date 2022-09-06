#ifndef CONTROLS_H
#define CONTROLS_H

#include <glm/glm.hpp>
#include <glm/detail/type_vec.hpp>

#if !defined(EMSCRIPTEN)
#include <glad/glad.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include "graphics.hpp"

void glfwScrollCallback(GLFWwindow* window, double xoffset, double yoffset);
#if defined(EMSCRIPTEN)
EM_BOOL emScrollCallback(int eventType, const EmscriptenWheelEvent *wheelEvent, void *userData);
#endif

void initControls(GLFWwindow* window);
void handleControls(GLFWwindow* window, Camera &camera, float dt);

namespace controls {
    extern glm::dvec2 scroll_offset;
    extern bool scrolled;
    extern bool left_mouse_click_release;
    extern bool left_mouse_click_press;

    extern glm::dvec2 mouse_position;
    extern glm::dvec2 delta_mouse_position;
}

#endif

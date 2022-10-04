#ifndef CONTROLS_H
#define CONTROLS_H

#include <glm/glm.hpp>

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include "graphics.hpp"

void glfwScrollCallback(GLFWwindow* window, double xoffset, double yoffset);
#if defined(EMSCRIPTEN)
EM_BOOL emScrollCallback(int eventType, const EmscriptenWheelEvent *wheelEvent, void *userData);
EM_BOOL emMouseDownCallback(int eventType, const EmscriptenMouseEvent *mouseEvent, void *userData);
EM_BOOL emMouseMoveCallback(int eventType, const EmscriptenMouseEvent *mouseEvent, void *userData);
EM_BOOL emMouseUpCallback(int eventType, const EmscriptenMouseEvent *mouseEvent, void *userData);
EM_BOOL emTouchStartCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData);
EM_BOOL emTouchMoveCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData);
EM_BOOL emTouchEndCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData);
EM_BOOL emTouchCancelCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData);
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

#ifndef GRAPHIC_HPP
#define GRAPHIC_HPP

#include <vector>

#if defined(EMSCRIPTEN)
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include <emscripten/html5.h>
#endif

#if !defined(EMSCRIPTEN)
#include <glad/glad.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include "assets.hpp"
#include "entities.hpp"

extern int    window_width;
extern int    window_height;
extern bool   window_resized;

struct Camera {
    const float near_plane = 0.5f, far_plane = 300.0f;
    const glm::vec3 up = glm::vec3(0,1,0);
    glm::vec3 position;
    glm::vec3 target;
    glm::mat4 view;
    glm::mat4 projection;
};

void createDefaultCamera(Camera &camera);
void updateCameraView(Camera &camera);
void updateCameraProjection(Camera &camera);
void updateCamera(Camera &camera);

#if defined(EMSCRIPTEN)
EM_BOOL emResizeCallback(int eventType, const EmscriptenUiEvent *uiEvent, void *userData);
#endif
void glfwResizeCallback(GLFWwindow* window, int width, int height);

void initGraphics();

void bindBackbuffer();
void drawEntities(const Entities &entities, const Camera &camera);

#endif

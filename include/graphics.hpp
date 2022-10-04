#ifndef GRAPHIC_HPP
#define GRAPHIC_HPP

#include <vector>

#if defined(EMSCRIPTEN)
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include <emscripten/html5.h>
#endif

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include "assets.hpp"
#include "entities.hpp"

struct Camera {
    const float near_plane = 0.5f, far_plane = 300.0f, fov = glm::radians(45.0f);
    bool is_ortho = false;

    const glm::vec3 up = glm::vec3(0,1,0);
    glm::vec3 position;
    glm::vec3 target;
    glm::mat4 view;
    glm::mat4 projection;

    // For now this selection is just the file were viewing
    bool selected = false;
    glm::vec3 selection_position;
    float selection_radius;  // Bounding sphere
};

void createDefaultCamera(Camera &camera);
void updateCameraView(Camera &camera);
void updateCameraProjection(Camera &camera);
void updateCamera(Camera &camera);
void positionCameraSelectionInView(Camera &camera, float fov_ratio);

#if !defined(EMSCRIPTEN)
void pushWriteFramebufferToTga(std::string path);
void checkWriteFrambufferToTga();
#endif

#if defined(EMSCRIPTEN)
EM_BOOL emResizeCallback(int eventType, const EmscriptenUiEvent *uiEvent, void *userData);
#endif
void glfwResizeCallback(GLFWwindow* window, int width, int height);

void initGraphics();

void bindBackbuffer();
void drawEntities(const Entities &entities, const Camera &camera);


extern int    window_width;
extern int    window_height;
extern bool   window_resized;

extern GLFWwindow *window;
namespace graphics{
    extern glm::vec3 sun_color;
    extern glm::vec3 sun_direction;
}

#endif

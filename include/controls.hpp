#ifndef CONTROLS_H
#define CONTROLS_H

#include <glm/glm.hpp>
#include "glm/detail/type_vec.hpp"

#include <GLFW/glfw3.h>

#include "globals.hpp"

void windowScrollCallback(GLFWwindow* window, double xoffset, double yoffset);

struct Camera {
    const float near_plane = 1.0f, far_plane = 100.0f;
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

namespace controls {
    extern glm::vec2 scroll_offset;
    extern bool scrolled;
    extern bool left_mouse_click_release;
    extern bool left_mouse_click_press;

    extern glm::dvec2 mouse_position;
    extern glm::dvec2 delta_mouse_position;
}
void handleControls(Camera &camera, float dt);
void initControls();

#endif

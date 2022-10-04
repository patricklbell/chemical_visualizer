#include <iostream>

#if !defined(EMSCRIPTEN)
#include "imgui.h"
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include "controls.hpp"
#include "graphics.hpp"
#include "utilities.hpp"

namespace controls {
    glm::dvec2 scroll_offset;
    bool scrolled = false;
    bool touched = false;

    bool left_mouse = false;
    bool right_mouse = false;
    bool mouse_moving = false;
    glm::dvec2 mouse_position;
    glm::dvec2 old_mouse_position;
}

using namespace controls;

void glfwScrollCallback(GLFWwindow* window, double xoffset, double yoffset){
    if(scroll_offset.x != xoffset || scroll_offset.y != yoffset) scrolled = true; 
    else                                                         scrolled = false;

    scroll_offset.x = xoffset;
    scroll_offset.y = yoffset;
}
#if defined(EMSCRIPTEN)
EM_BOOL emScrollCallback(int eventType, const EmscriptenWheelEvent *wheelEvent, void *userData)
{
    //std::cout << "Mousemove event, x: " << mouseEvent->clientX << ", y: " <<  mouseEvent->clientY  << "\n";
    double xoffset = wheelEvent->deltaX / 100.0f;
    double yoffset = wheelEvent->deltaY / 100.0f;
    if(scroll_offset.x != xoffset || scroll_offset.y != yoffset) scrolled = true; 
    else                                                         scrolled = false;

    scroll_offset.x = xoffset;
    scroll_offset.y = yoffset;
    return EM_TRUE;
}

EM_BOOL emMouseDownCallback(int eventType, const EmscriptenMouseEvent *mouseEvent, void *userData) {
    if(mouseEvent->button == 0)
        left_mouse = true;
    else if(mouseEvent->button == 2)
        right_mouse = true;
    return EM_TRUE;
}
EM_BOOL emMouseMoveCallback(int eventType, const EmscriptenMouseEvent *mouseEvent, void *userData) {
    old_mouse_position = mouse_position;
    mouse_position = glm::vec2(mouseEvent->targetX, mouseEvent->targetY);
    mouse_moving = true;
    return EM_TRUE;
}
EM_BOOL emMouseUpCallback(int eventType, const EmscriptenMouseEvent *mouseEvent, void *userData) {
    if(mouseEvent->button == 0)
        left_mouse = false;
    else if(mouseEvent->button == 2)
        right_mouse = false;
    return EM_TRUE;
}

EM_BOOL emTouchStartCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData) {
    touched = touchEvent->numTouches == 1;
    if(touched) {
        mouse_position = glm::vec2(touchEvent->touches[0].targetX, touchEvent->touches[0].targetY);
        old_mouse_position = mouse_position;
    }
    // Allow focus to return
    return EM_FALSE;
}
EM_BOOL emTouchMoveCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData) {
    if(touched) {
        old_mouse_position = mouse_position;
        mouse_position = glm::dvec2(touchEvent->touches[0].targetX, touchEvent->touches[0].targetY);
    }

    touched = touchEvent->numTouches == 1;
    return EM_FALSE;
}
EM_BOOL emTouchEndCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData) {
    touched = false;
    // Allow focus to return
    return EM_FALSE;
}
EM_BOOL emTouchCancelCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData) {
    touched = false;
    // Allow focus to return
    return EM_FALSE;
}

#endif

void initControls(GLFWwindow* window){
    glfwGetCursorPos(window, &mouse_position.x, &mouse_position.y);
    old_mouse_position = mouse_position;
}

void handleControls(GLFWwindow* window, Camera &camera, float dt) {
#if !defined(EMSCRIPTEN)
    auto old_mouse_position = mouse_position;
    glfwGetCursorPos(window, &mouse_position.x, &mouse_position.y);
    mouse_moving = mouse_position != old_mouse_position;
    left_mouse = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    right_mouse = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;

    ImGuiIO &io = ImGui::GetIO();
    if (io.WantCaptureMouse || io.WantCaptureKeyboard) return;
#endif

    auto delta_mouse_position = mouse_position - old_mouse_position;

    auto camera_look = camera.position - camera.target;
    auto camera_d = glm::length(camera_look);
    if ((left_mouse && mouse_moving) || touched) {
        mouse_moving = false;
        auto camera_right = glm::vec3(glm::transpose(camera.view)[0]);

        // Calculate the amount of rotation given the mouse movement.
        float delta_angle_x = (2 * glm::pi<float>() / (float)window_width); // a movement from left to right = 2PI
        float delta_angle_y = (glm::pi<float>() / (float)window_height);  // a movement from top to bottom = PI
        float x_angle = -delta_mouse_position.x * delta_angle_x;
        float y_angle = -delta_mouse_position.y * delta_angle_y;

        auto rotation_x = glm::angleAxis(x_angle, camera.up);
        camera_look = rotation_x * camera_look;

        // Handle camera passing over poles of orbit 
        // cos of angle between look and up is close to 1 -> parallel, -1 -> antiparallel
        auto l_cos_up = glm::dot(camera_look, camera.up) / camera_d;
        bool allow_rotation = true;
        if(abs(1 - l_cos_up) <= 0.01) {
            allow_rotation = y_angle > 0.f;
        } else if (abs(l_cos_up + 1) <= 0.01) {
            allow_rotation = y_angle < 0.f;
        }
        if (allow_rotation){
            auto rotation_y = glm::angleAxis(y_angle, camera_right);
            camera_look = rotation_y * camera_look;
        }

        // Update the camera view
        camera.position = camera_look + camera.target;
        updateCameraView(camera);
    } else if (right_mouse && mouse_moving) {
        auto inverse_vp = glm::inverse(camera.projection * camera.view);

            // Find Normalized Device coordinates mouse positions
            auto new_mouse_position_ndc = (mouse_position     / glm::dvec2((float)window_width, (float)window_height) - glm::dvec2(0.5)) * 2.0;
            auto old_mouse_position_ndc = (old_mouse_position / glm::dvec2((float)window_width, (float)window_height) - glm::dvec2(0.5)) * 2.0;

            // Project these mouse coordinates onto near plane to determine world coordinates
            auto new_mouse_position_world = glm::vec3(inverse_vp * glm::vec4(new_mouse_position_ndc.x, -new_mouse_position_ndc.y, 0, 1));
            auto old_mouse_position_world = glm::vec3(inverse_vp * glm::vec4(old_mouse_position_ndc.x, -old_mouse_position_ndc.y, 0, 1));

            // Scale movement such that point under mouse on plane of target (parallel to near plane) stays constant
            // Not needed for orthographic camera since it already exists in world space
            float ratio;
            if (camera.is_ortho) ratio = 1.0;
            else                 ratio = camera_d / camera.near_plane;
            auto delta = ratio * (new_mouse_position_world - old_mouse_position_world);

            camera.position -= delta;
            camera.target   -= delta; 
            updateCameraView(camera);
    }
    // Handles drift from last event when mouse button is held but not moving
    mouse_moving = false;

    if(scroll_offset.y != 0){
        float distance_scl = abs(1 + scroll_offset.y*0.1);
        if(distance_scl*camera_d < camera.far_plane - 10.0f) { // @hardcoded
            camera.position = camera.target + camera_look*distance_scl;
            updateCameraView(camera);
        }

        // Handle scroll event
        scroll_offset.y = 0;
    }
}

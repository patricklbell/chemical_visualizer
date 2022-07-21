#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// Include ImGui
#include "entities.hpp"
#include "glm/detail/func_geometric.hpp"
#include "glm/gtc/quaternion.hpp"
#include "imgui.h"

#include "controls.hpp"
#include "graphics.hpp"

namespace controls {
    glm::dvec2 scroll_offset;
    bool scrolled;
    glm::dvec2 mouse_position;
    glm::dvec2 delta_mouse_position;
}

using namespace controls;

void createDefaultCamera(Camera &camera){
    camera.position = glm::vec3(3,3,3);
    camera.target = glm::vec3(0,0,0);
    updateCameraView(camera);
    updateCameraProjection(camera);
}

void updateCameraView(Camera &camera){
    camera.view = glm::lookAt(camera.position, camera.target, camera.up);
}

void updateCameraProjection(Camera &camera){
    camera.projection = glm::perspective(glm::radians(45.0f), (float)window_width/(float)window_height, camera.near_plane, camera.far_plane);
}

void updateCamera(Camera &camera){
    updateCameraProjection(camera);
    updateCameraView(camera);
}

void windowScrollCallback(GLFWwindow* window, double xoffset, double yoffset){
    if(scroll_offset.x != xoffset || scroll_offset.y != yoffset) scrolled = true; 
    else                                                         scrolled = false;

    scroll_offset.x = xoffset;
    scroll_offset.y = yoffset;
}

void initControls(){
    glfwGetCursorPos(window, &mouse_position.x, &mouse_position.y);
    delta_mouse_position = glm::dvec2(0,0);
}

void handleControls(Camera &camera, float dt) {
    ImGuiIO& io = ImGui::GetIO();
    bool active = !(io.WantCaptureMouse || io.WantCaptureKeyboard);

    if(active){
        // Unlike other inputs, calculate delta but update mouse position immediately
        glm::dvec2 delta_mouse_position = mouse_position;
        glfwGetCursorPos(window, &mouse_position.x, &mouse_position.y);
        delta_mouse_position = mouse_position - delta_mouse_position;

        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            auto camera_right = glm::vec3(glm::transpose(camera.view)[0]);

            // Calculate the amount of rotation given the mouse movement.
            float delta_angle_x = (2 * PI / (float)window_width); // a movement from left to right = 2PI
            float delta_angle_y = (PI / (float)window_height);  // a movement from top to bottom = PI
            float x_angle = -delta_mouse_position.x * delta_angle_x;
            float y_angle = -delta_mouse_position.y * delta_angle_y;

            auto camera_look = camera.position - camera.target;

            auto rotation_x = glm::angleAxis(x_angle, camera.up);
            camera_look = rotation_x * camera_look;

            // Handle camera passing over poles of orbit 
            // cos of angle between look and up is close to 1 -> parallel, -1 -> antiparallel
            auto l_cos_up = glm::dot(camera_look, camera.up) / glm::length(camera_look);
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
        } else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
            auto camera_right = glm::vec3(glm::transpose(camera.view)[0]);

            auto delta = (float)delta_mouse_position.x*camera_right - (float)delta_mouse_position.y*camera.up;
            auto d = glm::length(camera.position - camera.target);
            camera.position -= 0.003f*d*delta;
            camera.target   -= 0.003f*d*delta;

            updateCameraView(camera);
        }

        if(scroll_offset.y != 0){
            float distance_scl = abs(1 + scroll_offset.y*0.1);

            camera.position = camera.target + (camera.position - camera.target)*distance_scl;
            updateCameraView(camera);

            // Handle scroll event
            scroll_offset.y = 0;
        }
    }    
}

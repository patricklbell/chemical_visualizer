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
#include "utilities.hpp"

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
    if (camera.is_ortho)
        updateCameraProjection(camera);
}

void positionCameraSelectionInView(Camera& camera, float fov_ratio) {
    camera.target = camera.selection_position;
    float min_screen_ratio = glm::min((float)window_height / (float)window_width, (float)window_width / (float)window_height);
    float d = camera.selection_radius / glm::max(glm::tan(min_screen_ratio * camera.fov * fov_ratio), 0.1f);
    d = glm::min(camera.far_plane / 2.0f, d);
    camera.position = camera.target + glm::normalize(-graphics::sun_direction) * d;
    updateCameraView(camera);
}

void updateCameraProjection(Camera &camera){
    if (camera.is_ortho) {
        auto distance = glm::length(camera.target - camera.position);
        auto ratio_size_per_depth = glm::atan(camera.fov);
        auto aspect = (float)window_width / (float)window_height;
        auto size_y = ratio_size_per_depth * distance;
        auto size_x = ratio_size_per_depth * aspect * distance;

        camera.projection = glm::ortho(-size_x, size_x, -size_y, size_y, -camera.near_plane, camera.far_plane);
    }
    else {
        camera.projection = glm::perspective(camera.fov, (float)window_width/(float)window_height, camera.near_plane, camera.far_plane);
    }
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
        glm::dvec2 old_mouse_position = mouse_position;
        glfwGetCursorPos(window, &mouse_position.x, &mouse_position.y);
        delta_mouse_position = (mouse_position - old_mouse_position);
        auto delta_mouse_position_n = delta_mouse_position / glm::dvec2((float)window_width, (float)window_height);

        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            auto camera_right = glm::vec3(glm::transpose(camera.view)[0]);

            // Calculate the amount of rotation given the mouse movement.
            float x_angle = -delta_mouse_position_n.x * 2 * PI; // a movement from left to right = 2PI
            float y_angle = -delta_mouse_position_n.y * PI; // a movement from top to bottom = PI

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
            else                 ratio = glm::length(camera.position - camera.target) / camera.near_plane;
            auto delta = ratio * (new_mouse_position_world - old_mouse_position_world);

            camera.position -= delta;
            camera.target   -= delta; 
            updateCameraView(camera);
        }

        if(scroll_offset.y != 0){
            float distance_scl = abs(1 + scroll_offset.y*0.1);

            camera.position = camera.target + (camera.position - camera.target) * distance_scl;

            updateCameraView(camera);

            // Handle scroll event
            scroll_offset.y = 0;
        }
    }    
}

#include <controls.hpp>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

#include <graphics.hpp>
#include <camera/core.hpp>

namespace Controls {
    glm::dvec2 scroll_offset;
    bool scrolled = false;
    bool touched  = false;

    bool left_mouse               = false;
    bool right_mouse              = false;
    bool mouse_moving             = false;
    glm::dvec2 mouse_position     = glm::vec2(-1, -1);
    glm::dvec2 old_mouse_position = glm::vec2(-1, -1);
}    // namespace Controls

using namespace Controls;

void handle_controls(Camera& camera, float dt) {
    auto delta_mouse_position = mouse_position - old_mouse_position;

    auto camera_d = glm::length(camera.forward);
    if ((left_mouse && mouse_moving) || touched) {
        mouse_moving = false;

        // Calculate the amount of rotation given the mouse movement.
        float x_screen_to_angle = (2 * glm::pi<float>() / (float)window_width);    // a movement from left to right = 2PI
        float y_screen_to_angle = (glm::pi<float>() / (float)window_height);       // a movement from top to bottom = PI
        float x_angle           = -delta_mouse_position.x * x_screen_to_angle;
        float y_angle           = -delta_mouse_position.y * y_screen_to_angle;

        auto rotation_x      = glm::angleAxis(x_angle, camera.up);
        auto rotated_forward = rotation_x * camera.forward;

        // Handle camera passing over poles of orbit
        // cos of angle between look and up is close to 1 -> parallel, -1 -> antiparallel
        auto l_cos_up       = glm::dot(rotated_forward, camera.up) / camera_d;
        bool allow_rotation = true;
        if (abs(1 - l_cos_up) <= 0.01) {
            allow_rotation = y_angle > 0.f;
        } else if (abs(l_cos_up + 1) <= 0.01) {
            allow_rotation = y_angle < 0.f;
        }
        if (allow_rotation) {
            auto rotation_y = glm::angleAxis(y_angle, camera.right);
            rotated_forward = rotation_y * rotated_forward;
        }

        camera.set_position(rotated_forward + camera.target);
    } else if (right_mouse && mouse_moving) {
        // Find Normalized Device coordinates mouse positions
        auto new_mouse_position_ndc = (mouse_position / glm::dvec2((float)window_width, (float)window_height) - glm::dvec2(0.5)) * 2.0;
        auto old_mouse_position_ndc = (old_mouse_position / glm::dvec2((float)window_width, (float)window_height) - glm::dvec2(0.5)) * 2.0;

        // Project these mouse coordinates onto near plane to determine world coordinates
        auto new_mouse_position_world = glm::vec3(camera.inv_vp * glm::vec4(new_mouse_position_ndc.x, -new_mouse_position_ndc.y, 0, 1));
        auto old_mouse_position_world = glm::vec3(camera.inv_vp * glm::vec4(old_mouse_position_ndc.x, -old_mouse_position_ndc.y, 0, 1));

        // Scale movement such that point under mouse on plane of target (parallel to near plane) stays constant
        // Not needed for orthographic camera since it already exists in world space
        float ratio = camera.type != Camera::Type::ORTHOGRAPHIC ? camera_d : 1.0;
        auto delta  = ratio * (new_mouse_position_world - old_mouse_position_world);

        camera.set_position(camera.position - delta);
        camera.set_target(camera.target - delta);
    }
    // Handles drift from last event when mouse button is held but not moving
    mouse_moving = false;

    if (scroll_offset.y != 0) {
        float distance_scl = abs(1 + scroll_offset.y * 0.1);
        distance_scl       = glm::clamp(distance_scl, 0.0f, camera.max_distance / camera_d);
        camera.set_position(camera.target + camera.forward * distance_scl);

        // Handle scroll event
        scroll_offset.y = 0;
    }
}

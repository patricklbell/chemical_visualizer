#ifndef CONTROLS_HPP
#define CONTROLS_HPP

#include <glm/glm.hpp>

#include <camera/core.hpp>

void handle_controls(Camera& camera, float dt);

namespace Controls {
    extern glm::dvec2 scroll_offset;
    extern bool scrolled;
    extern bool touched;

    extern bool left_mouse;
    extern bool right_mouse;
    extern bool mouse_moving;
    extern glm::dvec2 mouse_position;
    extern glm::dvec2 old_mouse_position;
}    // namespace Controls

#endif // CONTROLS_HPP

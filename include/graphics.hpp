#ifndef GRAPHICS_HPP
#define GRAPHICS_HPP

#include <camera/core.hpp>
#include <entities.hpp>

extern int window_width;
extern int window_height;

struct DrawSettings {
    glm::vec4 clear_color     = glm::vec4(0, 0, 0, 1);
    glm::vec3 light_direction = glm::vec3(0, -0.7071067811865475, -0.7071067811865475);
    glm::vec3 light_color     = 5.0f * glm::vec3(0.941, 0.933, 0.849);

    enum Mode {
        NORMAL = 0,
        GOOCH,
        BLINN_PHONG,
    } mode = Mode::NORMAL;
};

bool init_graphics();
void draw_entities(const Entities &entities, const Camera &camera, DrawSettings &draw_settings);

#ifndef EMSCRIPTEN
void screenshot_tga(std::string path);
#endif

#endif    // GRAPHICS_HPP
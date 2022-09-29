#ifndef GRAPHIC_HPP
#define GRAPHIC_HPP

#include <vector>

#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "controls.hpp"
#include "entities.hpp"

struct Entity;

extern int    window_width;
extern int    window_height;
extern bool   window_resized;

void windowSizeCallback(GLFWwindow* window, int width, int height);
void framebufferSizeCallback(GLFWwindow *window, int width, int height);

void initGraphicsPrimitives();

void bindBackbuffer();
void drawEntities(const Entities &entities, const Camera &camera);

void pushWriteFramebufferToTga(std::string_view path);
void checkWriteFrambufferToTga();

namespace graphics{
    extern Mesh sphere;
    extern Mesh cylinder;
}   
#endif

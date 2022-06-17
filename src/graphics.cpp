#include <algorithm>
#include <limits>
#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "glm/detail/type_mat.hpp"

#include "graphics.hpp"
#include "assets.hpp"
#include "entities.hpp"
#include "shader.hpp"
#include "utilities.hpp"
#include "shader.hpp"

int    window_width;
int    window_height;
bool   window_resized;
namespace graphics{
    glm::vec3 sun_color;
    glm::vec3 sun_direction;
    Mesh sphere;
    Mesh cylinder;
}

using namespace graphics;

void windowSizeCallback(GLFWwindow* window, int width, int height){
    if(width != window_width || height != window_height) window_resized = true;

    window_width  = width;
    window_height = height;

    // @note may not be necessary
    glViewport(0, 0, width, height);
}
void framebufferSizeCallback(GLFWwindow *window, int width, int height){}

void initGraphicsPrimitives(){
    sun_direction = glm::vec3(-0.7071067811865475, -0.7071067811865475, 0);
    sun_color = 10.0f*glm::vec3(0.941, 0.933, 0.849);

    loadMeshWithAssimp(sphere,   "data/models/sphere.obj");
    loadMeshWithAssimp(cylinder, "data/models/cylinder.obj");
}

void drawEntities(const std::vector<Entity*> &entities, const Camera &camera){
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    glClearColor(0.05,0.05,0.05,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(shader::basic_program);
    glUniform3fv(shader::basic_uniforms.sun_color, 1, &sun_color[0]);
    glUniform3fv(shader::basic_uniforms.sun_direction, 1, &sun_direction[0]);
    glUniform3fv(shader::basic_uniforms.camera_position, 1, &camera.position[0]);
    
    auto vp = camera.projection * camera.view;
    for (const auto &e : entities) {
        if(e->type != MESH_ENTITY) continue;
        auto m_e = (MeshEntity*)e;

        auto model = createModelMatrix(m_e->position, m_e->rotation, m_e->scale);
        auto mvp = vp * model;
        glUniformMatrix4fv(shader::basic_uniforms.mvp, 1, GL_FALSE, &mvp[0][0]);
        glUniformMatrix4fv(shader::basic_uniforms.model, 1, GL_FALSE, &model[0][0]);

        glUniform3fv(shader::basic_uniforms.albedo, 1, &m_e->albedo[0]);

        const auto &mesh = m_e->mesh;
        for (int j = 0; j < mesh->num_materials; ++j) {
            glBindVertexArray(mesh->vao);
            glDrawElements(mesh->draw_mode, mesh->draw_count[j], mesh->draw_type, (GLvoid*)(sizeof(GLubyte)*mesh->draw_start[j]));
        }
    }
}

void bindBackbuffer(){
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

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

    const Mesh *instanced_meshes[2] = {&graphics::cylinder, &graphics::sphere};
    std::vector<glm::mat4> instanced_models[2];
    std::vector<glm::mat4> instanced_mvps[2];
    std::vector<glm::vec3> instanced_albedos[2];
    
    auto vp = camera.projection * camera.view;
    for (const auto &e : entities) {
        if(e->type != MESH_ENTITY) continue;
        auto m_e = (MeshEntity*)e;

        // @speed precalculate since most entities are static
        auto model = createModelMatrix(m_e->position, m_e->rotation, m_e->scale);
        auto mvp = vp * model;
        
        for(int i = 0; i < 2; ++i){
            if(m_e->mesh == instanced_meshes[i]) {
                instanced_models[i].push_back(model);
                instanced_mvps[i].push_back(mvp);
                instanced_albedos[i].push_back(m_e->albedo);
            }
        }

        glUniformMatrix4fv(shader::basic_uniforms.mvp, 1, GL_FALSE, &mvp[0][0]);
        glUniformMatrix4fv(shader::basic_uniforms.model, 1, GL_FALSE, &model[0][0]);

        glUniform3fv(shader::basic_uniforms.albedo, 1, &m_e->albedo[0]);

        const auto &mesh = m_e->mesh;
        for (int j = 0; j < mesh->num_materials; ++j) {
            glBindVertexArray(mesh->vao);
            glDrawElements(mesh->draw_mode, mesh->draw_count[j], mesh->draw_type, (GLvoid*)(sizeof(GLubyte)*mesh->draw_start[j]));
        }
    }

    glUseProgram(shader::basic_instanced_program);
    glUniform3fv(shader::basic_instanced_uniforms.sun_color, 1, &sun_color[0]);
    glUniform3fv(shader::basic_instanced_uniforms.sun_direction, 1, &sun_direction[0]);
    glUniform3fv(shader::basic_instanced_uniforms.camera_position, 1, &camera.position[0]);

    // @speed I think this could be done by making one shared buffer for all instanced meshes
    // and indexing into it which would propably be faster
    for(int i = 0; i < 2; ++i) {
        auto &mesh = instanced_meshes[i];
        int num = instanced_albedos[i].size();
        glBindVertexArray(mesh->vao);

        unsigned int albedo_vbo;
        glGenBuffers(1, &albedo_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, albedo_vbo);
        glBufferData(GL_ARRAY_BUFFER, num * sizeof(glm::vec3), &instanced_albedos[i][0], GL_STATIC_DRAW);
        glVertexAttribPointer(2, 3, GL_FLOAT, false, 0, 0);
        glEnableVertexAttribArray(2);

        constexpr auto v4_s = sizeof(glm::vec4);

        unsigned int model_vbo;
        glGenBuffers(1, &model_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, model_vbo);
        glBufferData(GL_ARRAY_BUFFER, num*sizeof(glm::mat4), &instanced_models[i][0], GL_STATIC_DRAW);

        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4*v4_s, (void*)0);
        glEnableVertexAttribArray(4);
        glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4*v4_s, (void*)(1*v4_s));
        glEnableVertexAttribArray(5);
        glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4*v4_s, (void*)(2*v4_s));
        glEnableVertexAttribArray(6);
        glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, 4*v4_s, (void*)(3*v4_s));
        glVertexAttribDivisor(3, 1);
        glVertexAttribDivisor(4, 1);
        glVertexAttribDivisor(5, 1);
        glVertexAttribDivisor(6, 1);

        unsigned int mvp_vbo;
        glGenBuffers(1, &mvp_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, mvp_vbo);
        glBufferData(GL_ARRAY_BUFFER, num*sizeof(glm::mat4), &instanced_mvps[i][0], GL_STATIC_DRAW);

        glEnableVertexAttribArray(7);
        glVertexAttribPointer(7, 4,  GL_FLOAT, GL_FALSE, 4*v4_s, (void*)0);
        glEnableVertexAttribArray(8);
        glVertexAttribPointer(8, 4,  GL_FLOAT, GL_FALSE, 4*v4_s, (void*)(1*v4_s));
        glEnableVertexAttribArray(9);
        glVertexAttribPointer(9, 4,  GL_FLOAT, GL_FALSE, 4*v4_s, (void*)(2*v4_s));
        glEnableVertexAttribArray(10);
        glVertexAttribPointer(10, 4, GL_FLOAT, GL_FALSE, 4*v4_s, (void*)(3*v4_s));
        glVertexAttribDivisor(7,  1);
        glVertexAttribDivisor(8,  1);
        glVertexAttribDivisor(9,  1);
        glVertexAttribDivisor(10, 1);

        // @assumption instanced meshes have only one submesh/material
        glDrawElementsInstanced(mesh->draw_mode, mesh->draw_count[0], mesh->draw_type, (GLvoid*)(sizeof(GLubyte)*mesh->draw_start[0]), num);
    }
}

void bindBackbuffer(){
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

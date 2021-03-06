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
}

void drawEntities(const Entities &entities, const Camera &camera){
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    glClearColor(0.12,0.13, 0.2,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(shader::basic_program);
    glUniform3fv(shader::basic_uniforms.sun_color, 1, &sun_color[0]);
    glUniform3fv(shader::basic_uniforms.sun_direction, 1, &sun_direction[0]);
    glUniform3fv(shader::basic_uniforms.camera_position, 1, &camera.position[0]);

    auto vp = camera.projection * camera.view;
    for (const auto &m_e : entities.mesh_entities) {
        // @speed precalculate since most entities are static
        auto model = createModelMatrix(m_e.position, m_e.rotation, m_e.scale);
        auto mvp = vp * model;
        
        glUniformMatrix4fv(shader::basic_uniforms.mvp, 1, GL_FALSE, &mvp[0][0]);
        glUniformMatrix4fv(shader::basic_uniforms.model, 1, GL_FALSE, &model[0][0]);

        glUniform3fv(shader::basic_uniforms.albedo, 1, &m_e.albedo[0]);

        const auto &mesh = m_e.mesh;
        glBindVertexArray(mesh.vao);
        for (int j = 0; j < mesh.num_materials; ++j) {
            glDrawElements(mesh.draw_mode, mesh.draw_count[j], mesh.draw_type, (GLvoid*)(sizeof(GLubyte)*mesh.draw_start[j]));
        }
    }

    int num_instance_meshes = entities.instanced_meshes.size();
    int num_instance_entities = entities.instanced_entities.size();
    if(num_instance_meshes > 0 && num_instance_entities > 0) {
        std::vector<std::vector<glm::mat4>> instanced_models(num_instance_meshes, std::vector<glm::mat4>());
        std::vector<std::vector<glm::mat4>> instanced_mvps(num_instance_meshes, std::vector<glm::mat4>());
        std::vector<std::vector<glm::vec3>> instanced_albedos(num_instance_meshes, std::vector<glm::vec3>());
    
        for(const auto &i_e : entities.instanced_entities) {
            auto model = createModelMatrix(i_e.position, i_e.rotation, i_e.scale);
            auto mvp = vp * model;

            instanced_models[i_e.instance_mesh].push_back(model);
            instanced_mvps[i_e.instance_mesh].push_back(mvp);
            instanced_albedos[i_e.instance_mesh].push_back(i_e.albedo);
        }

        glUseProgram(shader::basic_instanced_program);
        glUniform3fv(shader::basic_instanced_uniforms.sun_color, 1, &sun_color[0]);
        glUniform3fv(shader::basic_instanced_uniforms.sun_direction, 1, &sun_direction[0]);
        glUniform3fv(shader::basic_instanced_uniforms.camera_position, 1, &camera.position[0]);

        // @speed I think this could be done by making one shared buffer for all instanced meshes
        // and indexing into it which would propably be faster
        for(int i = 0; i < num_instance_meshes; ++i) {
            auto &mesh = entities.instanced_meshes[i];

            int num = instanced_albedos[i].size();
            glBindVertexArray(mesh.vao);

            GLuint albedo_vbo;
            glGenBuffers(1, &albedo_vbo);
            glBindBuffer(GL_ARRAY_BUFFER, albedo_vbo);
            glBufferData(GL_ARRAY_BUFFER, num*sizeof(glm::vec3), &instanced_albedos[i][0], GL_STATIC_DRAW);

            glEnableVertexAttribArray(2);
            glVertexAttribPointer(2, 3, GL_FLOAT, false, 0, 0);
            glVertexAttribDivisor(2, 1);  

            GLuint model_vbo;
            glGenBuffers(1, &model_vbo);
            glBindBuffer(GL_ARRAY_BUFFER, model_vbo);
            glBufferData(GL_ARRAY_BUFFER, num*sizeof(glm::mat4), &instanced_models[i][0], GL_STATIC_DRAW);

            constexpr auto v4_s = sizeof(glm::vec4);
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

            GLuint mvp_vbo;
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

            glBindVertexArray(mesh.vao);
            for (int j = 0; j < mesh.num_materials; ++j) {
                glDrawElementsInstanced(mesh.draw_mode, mesh.draw_count[j], mesh.draw_type, (GLvoid*)(sizeof(GLubyte)*mesh.draw_start[j]), num);
            }
        }
    }
}

void bindBackbuffer(){
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

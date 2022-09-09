#include <algorithm>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#if !defined(EMSCRIPTEN)
#include <glad/glad.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "graphics.hpp"
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
    camera.projection = glm::perspective(camera.fov, (float)window_width/(float)window_height, camera.near_plane, camera.far_plane);
}

void updateCamera(Camera &camera){
    updateCameraProjection(camera);
    updateCameraView(camera);
}

using namespace graphics;

#if defined(EMSCRIPTEN)
EM_BOOL emResizeCallback(int eventType, const EmscriptenUiEvent *uiEvent, void *userData)
{
    int width = uiEvent->windowInnerWidth;
    int height = uiEvent->windowInnerHeight;
#if DEBUG
    std::cout << "Emscripten resize event, width: " << width << ", height: " <<  height << "\n";
#endif

    if(width != window_width || height != window_height) window_resized = true;
    window_width  = width;
    window_height = height;

    // @note may not be necessary
    glViewport(0, 0, width, height);

    return EM_TRUE;
}
#endif

void glfwResizeCallback(GLFWwindow* window, int width, int height) {
    std::cout << "Glfw resize event, width: " << width << ", height: " <<  height << "\n";

    if(width != window_width || height != window_height) window_resized = true;

    window_width  = width;
    window_height = height;

    // @note may not be necessary
    glViewport(0, 0, width, height);
}

void initGraphics(){
    sun_direction = glm::vec3(-0.7071067811865475, -0.7071067811865475, 0);
    sun_color = 5.0f*glm::vec3(0.941, 0.933, 0.849);

    readMeshFile(sphere, "data/models/sphere.mesh");
    readMeshFile(cylinder, "data/models/cylinder.mesh");
}

void drawEntities(const Entities &entities, const Camera &camera){
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    glClearColor(0.12,0.13, 0.2, 1);
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

    int num_instance_entities = entities.instanced_entities.size();
    constexpr int num_instance_meshes = (int)InstanceNames::NUM;
    if(num_instance_entities > 0) {
        std::array<std::vector<glm::mat4>, num_instance_meshes> instanced_models{std::vector<glm::mat4>()};
        std::array<std::vector<glm::mat4>, num_instance_meshes> instanced_mvps{std::vector<glm::mat4>()};
        std::array<std::vector<glm::vec3>, num_instance_meshes> instanced_albedos{std::vector<glm::vec3>()};
    
        for(const auto &i_e : entities.instanced_entities) {
            const int i = (int)i_e.instance_mesh;
            if(i < 0 || i >= num_instance_meshes) continue;

            auto model = createModelMatrix(i_e.position, i_e.rotation, i_e.scale);
            auto mvp = vp * model;

            instanced_models[i].push_back(model);
            instanced_mvps[i].push_back(mvp);
            instanced_albedos[i].push_back(i_e.albedo);
        }

        glUseProgram(shader::basic_instanced_program);
        glUniform3fv(shader::basic_instanced_uniforms.sun_color, 1, &sun_color[0]);
        glUniform3fv(shader::basic_instanced_uniforms.sun_direction, 1, &sun_direction[0]);
        glUniform3fv(shader::basic_instanced_uniforms.camera_position, 1, &camera.position[0]);

        // @speed I think this could be done by making one shared buffer for all instanced meshes
        // and indexing into it which would propably be faster
        for(int i = 0; i < num_instance_meshes; ++i) {
            Mesh *mesh;
            switch(InstanceNames(i)) {
                case InstanceNames::SPHERE:
                    mesh = &sphere;
                    break;
                case InstanceNames::CYLINDER:
                    mesh = &cylinder;
                    break;
                default:
                    continue;
            }

            int num = instanced_albedos[i].size();
            glBindVertexArray(mesh->vao);

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

            glBindVertexArray(mesh->vao);
            for (int j = 0; j < mesh->num_materials; ++j) {
                glDrawElementsInstanced(mesh->draw_mode, mesh->draw_count[j], mesh->draw_type, (GLvoid*)(sizeof(*mesh->indices)*mesh->draw_start[j]), num);
            }
        }
    }
}

void bindBackbuffer(){
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

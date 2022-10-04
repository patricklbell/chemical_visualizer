#include <algorithm>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#if defined(EMSCRIPTEN)
#include <emscripten.h>
#endif

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "glm/detail/type_mat.hpp"

#include "graphics.hpp"
#include "utilities.hpp"
#include "shader.hpp"
#include "ui.hpp"

int    window_width;
int    window_height;
bool   window_resized;
GLFWwindow *window;

namespace graphics{
    glm::vec3 sun_color;
    glm::vec3 sun_direction;
    
    std::vector<Mesh> instanced_meshes;
}
using namespace graphics;

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

    float aspect = (float)window_width / (float)window_height;
    float fov_x  = glm::atan(glm::tan(camera.fov / 2.0f) * aspect) * 2.0f;
    float min_fov = glm::min(fov_x, camera.fov);

    float d = camera.selection_radius / glm::tan(min_fov * fov_ratio / 2.0f);
    d = glm::min(camera.far_plane / 1.5f, d);
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

// @todo make screenshots work for emscripten
#if !defined(EMSCRIPTEN)
bool do_write_tga = false;
std::string tga_path = "";

void pushWriteFramebufferToTga(std::string path)
{
    do_write_tga = true;
    tga_path = path;
}

void checkWriteFrambufferToTga()
{
    if (!do_write_tga) return;

    int *buffer = new int[window_width * window_height * 3];
    glReadPixels(0, 0, window_width, window_height, GL_BGR, GL_UNSIGNED_BYTE,
                 buffer);
    FILE *fp = fopen(tga_path.data(), "wb");
    if (fp == NULL)
    {
        fprintf(stderr, "Failed to open file %s to write TGA.",
                tga_path.data());
        return;
    }

    printf("----------------Writing Frame to TGA %s----------------\n",
           tga_path.data());
    short TGAhead[] = {
        0, 2, 0, 0, 0, 0, (short)window_width, (short)window_height, 24};
    fwrite(TGAhead, sizeof(TGAhead), 1, fp);
    fwrite(buffer, 3 * window_width * window_height, 1, fp);
    delete[] buffer;

    fclose(fp);

    do_write_tga = false;
}
#endif

#if defined(EMSCRIPTEN)
EM_BOOL emResizeCallback(int eventType, const EmscriptenUiEvent *uiEvent, void *userData)
{
    int width, height;
    if(ui::fullscreen) {
        width = uiEvent->windowInnerWidth;
        height = uiEvent->windowInnerHeight;
    } else {
        double w,h;
        emscripten_get_element_css_size("#canvas", &w, &h);
        width = (int)w;
        height = (int)h; 
    }
#if DEBUG
    std::cout << "Emscripten resize event, width: " << width << ", height: " <<  height << "\n";
#endif

    if(width != window_width || height != window_height) window_resized = true;
    window_width  = width;
    window_height = height;

    glfwSetWindowSize(window, window_width, window_height);
    if(window_resized) {
        emscripten_set_element_css_size("#canvas", (double)window_width, (double)window_height);
    }

    return EM_TRUE;
}
#endif

void glfwResizeCallback(GLFWwindow* window, int width, int height) {
#if DEBUG
    std::cout << "Glfw resize event, width: " << width << ", height: " <<  height << "\n";
#endif

    if(width != window_width || height != window_height) window_resized = true;

    window_width  = width;
    window_height = height;

    // @note may not be necessary
    glViewport(0, 0, width, height);
}

void initGraphics(){
    sun_direction = glm::vec3(0, -0.7071067811865475, -0.7071067811865475);
    sun_color = 5.0f*glm::vec3(0.941, 0.933, 0.849);

    instanced_meshes.resize(InstancedMeshType::NUM_TYPES);
#ifdef USE_ASSIMP
    loadMeshWithAssimp(instanced_meshes[InstancedMeshType::SPHERE],
                       "data/models/sphere.obj");
    loadMeshWithAssimp(instanced_meshes[InstancedMeshType::CYLINDER],
                       "data/models/cylinder.obj");
#else
    readMeshFile(instanced_meshes[InstancedMeshType::SPHERE],
                 "data/models/sphere.mesh");
    readMeshFile(instanced_meshes[InstancedMeshType::CYLINDER],
                 "data/models/cylinder.mesh");
#endif
}

void drawEntities(const Entities &entities, const Camera &camera){
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    if (ui::dark_mode)
        glClearColor(0.12, 0.13, 0.2, 1.0);
    else
        glClearColor(0.95, 0.95, 0.95, 1.0);
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

        glUniform4fv(shader::basic_uniforms.albedo, 1, &m_e.albedo[0]);

        const auto mesh = m_e.mesh;
        glBindVertexArray(mesh->vao);
        for (int j = 0; j < mesh->num_materials; ++j) {
            glDrawElements(mesh->draw_mode, mesh->draw_count[j], mesh->draw_type, (GLvoid*)(sizeof(GLubyte)*mesh->draw_start[j]));
        }
    }

    int num_instance_meshes = instanced_meshes.size();
    int num_instance_entities = entities.instanced_entities.size();
    if(num_instance_entities > 0) {
        std::vector<std::vector<glm::mat4>> instanced_models(num_instance_meshes, std::vector<glm::mat4>());
        std::vector<std::vector<glm::mat4>> instanced_mvps(num_instance_meshes, std::vector<glm::mat4>());
        std::vector<std::vector<glm::vec4>> instanced_albedos(num_instance_meshes, std::vector<glm::vec4>());
    
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
            auto &mesh = instanced_meshes[i];

            int num = instanced_albedos[i].size();
            glBindVertexArray(mesh.vao);

            GLuint albedo_vbo;
            glGenBuffers(1, &albedo_vbo);
            glBindBuffer(GL_ARRAY_BUFFER, albedo_vbo);
            glBufferData(GL_ARRAY_BUFFER, num*sizeof(instanced_albedos[0][0]), &instanced_albedos[i][0], GL_STATIC_DRAW);

            glEnableVertexAttribArray(2);
            glVertexAttribPointer(2, 4, GL_FLOAT, false, 0, 0);
            glVertexAttribDivisor(2, 1);  

            GLuint model_vbo;
            glGenBuffers(1, &model_vbo);
            glBindBuffer(GL_ARRAY_BUFFER, model_vbo);
            glBufferData(GL_ARRAY_BUFFER, num*sizeof(instanced_models[0][0]), &instanced_models[i][0], GL_STATIC_DRAW);

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
            glBufferData(GL_ARRAY_BUFFER, num*sizeof(instanced_mvps[0][0]), &instanced_mvps[i][0], GL_STATIC_DRAW);

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
                glDrawElementsInstanced(mesh.draw_mode, mesh.draw_count[j], mesh.draw_type, (GLvoid*)(sizeof(*mesh.indices)*mesh.draw_start[j]), num);
            }

            glDeleteBuffers(1, &albedo_vbo);
            glDeleteBuffers(1, &model_vbo);
            glDeleteBuffers(1, &mvp_vbo);
        }
    }

    glDisable(GL_BLEND);
}

void bindBackbuffer(){
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

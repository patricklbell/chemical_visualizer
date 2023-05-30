#include <graphics.hpp>

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <limits>

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif

#include <glm/glm.hpp>

#include <shader/globals.hpp>
#include <camera/core.hpp>
#include <utilities/math.hpp>
#include <mesh.hpp>
#include <entities.hpp>

int window_width  = 100;
int window_height = 100;

#ifndef EMSCRIPTEN
void do_screenshots();
#endif

std::vector<Mesh> instanced_meshes;
std::vector<std::vector<glm::mat4>> instanced_models;
std::vector<std::vector<glm::mat4>> instanced_mvps;
std::vector<std::vector<glm::vec4>> instanced_albedos;

[[nodiscard]] bool init_graphics() {
    bool failed = !init_shaders();

    instanced_meshes.resize(InstancedEntity::Type::NUM);
    failed |= !instanced_meshes[InstancedEntity::Type::SPHERE].load("data/models/sphere.mesh");
    failed |= !instanced_meshes[InstancedEntity::Type::CYLINDER].load("data/models/cylinder.mesh");

    instanced_models  = std::vector<std::vector<glm::mat4>>(InstancedEntity::Type::NUM, std::vector<glm::mat4>());
    instanced_mvps    = std::vector<std::vector<glm::mat4>>(InstancedEntity::Type::NUM, std::vector<glm::mat4>());
    instanced_albedos = std::vector<std::vector<glm::vec4>>(InstancedEntity::Type::NUM, std::vector<glm::vec4>());

    return !failed;
}

void draw_entities(const Entities &entities, const Camera &camera, DrawSettings &draw_settings) {
    // These state calls are not really necessary, @todo proper gl state management
    glViewport(0, 0, window_width, window_height);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(draw_settings.clear_color.x, draw_settings.clear_color.y, draw_settings.clear_color.z, draw_settings.clear_color.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    Shaders::basic.set_macro("PBR", draw_settings.mode == DrawSettings::Mode::NORMAL, false);
    Shaders::basic.set_macro("GOOCH", draw_settings.mode == DrawSettings::Mode::GOOCH, false);
    Shaders::basic.set_macro("DIRECTIONAL", draw_settings.mode == DrawSettings::Mode::BLINN_PHONG, false);
    Shaders::basic.apply_macros();
    Shaders::basic_instanced.set_macro("PBR", draw_settings.mode == DrawSettings::Mode::NORMAL, false);
    Shaders::basic_instanced.set_macro("GOOCH", draw_settings.mode == DrawSettings::Mode::GOOCH, false);
    Shaders::basic_instanced.set_macro("DIRECTIONAL", draw_settings.mode == DrawSettings::Mode::BLINN_PHONG, false);
    Shaders::basic_instanced.apply_macros();

    Shaders::basic.bind();
    Shaders::basic.uniform("light_color", draw_settings.light_color);
    Shaders::basic.uniform("light_direction", draw_settings.light_direction);
    Shaders::basic.uniform("camera_position", camera.position);

    for (const auto &m_e : entities.mesh_entities) {
        // @speed precalculate since most entities are static
        auto model = model_matrix(m_e.position, m_e.rotation, m_e.scale);
        auto mvp   = camera.vp * model;

        Shaders::basic.uniform("mvp", mvp);
        Shaders::basic.uniform("model", model);
        Shaders::basic.uniform("albedo", m_e.albedo);

        const auto mesh = m_e.mesh;
        glBindVertexArray(mesh->vao);
        for (int j = 0; j < mesh->num_submeshes; ++j) {
            glDrawElements(mesh->draw_mode, mesh->draw_count[j],
                           mesh->draw_type,
                           (GLvoid *)(sizeof(GLubyte) * mesh->draw_start[j]));
        }
    }

    int num_instance_meshes   = instanced_meshes.size();
    int num_instance_entities = entities.instanced_entities.size();
    if (num_instance_entities > 0) {
        for (int i = 0; i < InstancedEntity::Type::NUM; i++) {
            instanced_models[i].clear();
            instanced_mvps[i].clear();
            instanced_albedos[i].clear();
        }

        for (const auto &i_e : entities.instanced_entities) {
            const int i = (int)i_e.type;
            if (i < 0 || i >= num_instance_meshes)
                continue;

            auto model = model_matrix(i_e.position, i_e.rotation, i_e.scale);
            auto mvp   = camera.vp * model;

            instanced_models[i].push_back(model);
            instanced_mvps[i].push_back(mvp);
            instanced_albedos[i].push_back(i_e.albedo);
        }

        Shaders::basic_instanced.bind();
        Shaders::basic_instanced.uniform("light_color", draw_settings.light_color);
        Shaders::basic_instanced.uniform("light_direction", draw_settings.light_direction);
        Shaders::basic_instanced.uniform("camera_position", camera.position);

        // @speed I think this could be done by making one shared buffer for all
        // instanced meshes and indexing into it which would propably be faster
        for (int i = 0; i < num_instance_meshes; ++i) {
            auto &mesh = instanced_meshes[i];

            int num_submeshes = instanced_albedos[i].size();
            glBindVertexArray(mesh.vao);

            GLuint albedo_vbo;
            glGenBuffers(1, &albedo_vbo);
            glBindBuffer(GL_ARRAY_BUFFER, albedo_vbo);
            glBufferData(GL_ARRAY_BUFFER, num_submeshes * sizeof(instanced_albedos[0][0]),
                         &instanced_albedos[i][0], GL_STATIC_DRAW);

            glEnableVertexAttribArray(2);
            glVertexAttribPointer(2, 4, GL_FLOAT, false, 0, 0);
            glVertexAttribDivisor(2, 1);

            GLuint model_vbo;
            glGenBuffers(1, &model_vbo);
            glBindBuffer(GL_ARRAY_BUFFER, model_vbo);
            glBufferData(GL_ARRAY_BUFFER, num_submeshes * sizeof(instanced_models[0][0]),
                         &instanced_models[i][0], GL_STATIC_DRAW);

            constexpr auto v4_s = sizeof(glm::vec4);
            glEnableVertexAttribArray(3);
            glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * v4_s,
                                  (void *)0);
            glEnableVertexAttribArray(4);
            glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 4 * v4_s,
                                  (void *)(1 * v4_s));
            glEnableVertexAttribArray(5);
            glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 4 * v4_s,
                                  (void *)(2 * v4_s));
            glEnableVertexAttribArray(6);
            glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, 4 * v4_s,
                                  (void *)(3 * v4_s));
            glVertexAttribDivisor(3, 1);
            glVertexAttribDivisor(4, 1);
            glVertexAttribDivisor(5, 1);
            glVertexAttribDivisor(6, 1);

            GLuint mvp_vbo;
            glGenBuffers(1, &mvp_vbo);
            glBindBuffer(GL_ARRAY_BUFFER, mvp_vbo);
            glBufferData(GL_ARRAY_BUFFER, num_submeshes * sizeof(instanced_mvps[0][0]),
                         &instanced_mvps[i][0], GL_STATIC_DRAW);

            glEnableVertexAttribArray(7);
            glVertexAttribPointer(7, 4, GL_FLOAT, GL_FALSE, 4 * v4_s,
                                  (void *)0);
            glEnableVertexAttribArray(8);
            glVertexAttribPointer(8, 4, GL_FLOAT, GL_FALSE, 4 * v4_s,
                                  (void *)(1 * v4_s));
            glEnableVertexAttribArray(9);
            glVertexAttribPointer(9, 4, GL_FLOAT, GL_FALSE, 4 * v4_s,
                                  (void *)(2 * v4_s));
            glEnableVertexAttribArray(10);
            glVertexAttribPointer(10, 4, GL_FLOAT, GL_FALSE, 4 * v4_s,
                                  (void *)(3 * v4_s));
            glVertexAttribDivisor(7, 1);
            glVertexAttribDivisor(8, 1);
            glVertexAttribDivisor(9, 1);
            glVertexAttribDivisor(10, 1);

            glBindVertexArray(mesh.vao);
            for (int j = 0; j < mesh.num_submeshes; ++j) {
                glDrawElementsInstanced(
                    mesh.draw_mode, mesh.draw_count[j], mesh.draw_type,
                    (GLvoid *)(sizeof(*mesh.indices) * mesh.draw_start[j]),
                    num_submeshes);
            }

            glDeleteBuffers(1, &albedo_vbo);
            glDeleteBuffers(1, &model_vbo);
            glDeleteBuffers(1, &mvp_vbo);
        }
    }

#ifndef EMSCRIPTEN
    do_screenshots();
#endif
}

#ifndef EMSCRIPTEN
std::string tga_path = "";

void screenshot_tga(std::string_view path) {
    tga_path = path;
}

void do_screenshots() {
    if (tga_path == "")
        return;

    int *buffer = new int[window_width * window_height * 3];
    glReadPixels(0, 0, window_width, window_height, GL_BGR, GL_UNSIGNED_BYTE,
                 buffer);
    FILE *fp = fopen(tga_path.data(), "wb");
    if (fp == NULL) {
        std::cerr << "Failed to open file " << tga_path << " to write TGA.\n";
        return;
    }

    short TGAhead[] = {
        0, 2, 0, 0, 0, 0, (short)window_width, (short)window_height, 24};
    fwrite(TGAhead, sizeof(TGAhead), 1, fp);
    fwrite(buffer, 3 * window_width * window_height, 1, fp);
    delete[] buffer;

    fclose(fp);

#ifdef VERBOSE
    std::cout << "Wrote framebuffer as TGA to " << tga_path << "\n";
#endif

    tga_path = "";
}

#endif    // !EMSCRIPTEN

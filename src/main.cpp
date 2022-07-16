#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <stack>
#include <array>
#include <map>
#include <filesystem>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>
#ifdef _WINDOWS
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3native.h>
#include <ShellScalingApi.h>
#endif
GLFWwindow* window;

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "globals.hpp"
#include "shader.hpp"
#include "texture.hpp"
#include "utilities.hpp"
#include "graphics.hpp"
#include "controls.hpp"
#include "editor.hpp"
#include "assets.hpp"
#include "entities.hpp"
#include "loader.hpp"

int main() {
    if(!glfwInit()) {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return 0;
    }
    
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    // GL 4.1 + GLSL 410
    shader::glsl_version = "#version 410\n";
    glfwWindowHint( // required on Mac OS
        GLFW_OPENGL_FORWARD_COMPAT,
        GL_TRUE
    );
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
#elif __linux__
    // GL 4.3 + GLSL 430
    shader::glsl_version = "#version 430\n";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#elif _WIN32
    // GL 3.0 + GLSL 130
    shader::glsl_version = "#version 130\n";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
#else
    shader::glsl_version = "#version 130\n";
#endif

#if __APPLE__
    // to prevent 1200x800 from becoming 2400x1600
    glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_FALSE);
#endif

    window = glfwCreateWindow(1024, 700, "Window", NULL, NULL);
    if(window == NULL) {
        fprintf( stderr, "Failed to open GLFW window.\n" );
        getchar();
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);

    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        glfwTerminate();
        return 1;
    }

    // Ensure we can capture keys correctly 
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

    glfwSetScrollCallback(window, windowScrollCallback);
    glfwGetWindowSize(window, &window_width, &window_height);
    glfwSetWindowSizeCallback(window, windowSizeCallback);

    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
    glfwSwapInterval(0);

    glEnable(GL_MULTISAMPLE);

    initGraphicsPrimitives();

    Camera camera;
    createDefaultCamera(camera);

    // Load shaders
    loadBasicShader("data/shaders/basic.gl");
    loadBasicInstancedShader("data/shaders/basic_instanced.gl");

    // Map for last update time for hotswaping files
    std::filesystem::file_time_type empty_file_time;
    std::map<const std::string, std::pair<shader::TYPE, std::filesystem::file_time_type>> shader_update_times = {
        {"data/shaders/basic.gl", {shader::TYPE::BASIC_SHADER, empty_file_time}},
        {"data/shaders/basic_instanced.gl", {shader::TYPE::BASIC_INSTANCED_SHADER, empty_file_time}},
    };
    // Fill in with correct file time
    for (auto &pair : shader_update_times) {
        if(std::filesystem::exists(pair.first)) 
            pair.second.second = std::filesystem::last_write_time(pair.first);
    }

    PdbFile pdb_file;
    loadPdbFile(pdb_file, "data/examples/pdb/1bzv.pdb");
   
    Entities entities;
    camera.target = createEntitiesFromPdbFile(entities, pdb_file);

    initGui();
    initControls();

    double last_time = glfwGetTime();
    double last_filesystem_hotswap_check = last_time;
    window_resized = true;
    do {
        double current_time = glfwGetTime();
        float true_dt = current_time - last_time;
        last_time = current_time;
        static const float dt = 1.0/60.0;

        if (window_resized){
            updateCamera(camera);
        }

        // @developer Hotswap shader files
        if(current_time - last_filesystem_hotswap_check >= 2.0){
            last_filesystem_hotswap_check = current_time;
            for (auto &pair : shader_update_times) {
                if(pair.second.second != std::filesystem::last_write_time(pair.first)){
                    pair.second.second = std::filesystem::last_write_time(pair.first);
                    switch (pair.second.first) {
                        case shader::TYPE::BASIC_SHADER:
                            loadBasicShader(pair.first.c_str());
                            break;
                        case shader::TYPE::BASIC_INSTANCED_SHADER:
                            loadBasicInstancedShader(pair.first.c_str());
                            break;
                        default:
                            break;
                    }
                } 
            }
        }
        handleControls(camera, true_dt);

        bindBackbuffer();
        drawEntities(entities, camera);

        drawGui(camera, entities);

        // swap backbuffer with front buffer
        glfwSwapBuffers(window);
        window_resized = false;
        glfwPollEvents();

        GLenum code = glGetError();
        if(code != GL_NO_ERROR){
            static const GLubyte* err;
            err = gluErrorString(code);
            fprintf(stderr, "<--------------------OpenGL ERROR-------------------->\n%s\n", err);
        }
    } // exit if the ESC key was pressed or the window was closed
    while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS && glfwWindowShouldClose(window) == 0 );

    deleteShaderPrograms();    

    glfwTerminate();
    return 0;
}



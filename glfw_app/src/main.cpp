#include <array>
#include <chrono>
#include <iostream>
#include <map>
#include <thread>

#if defined(EMSCRIPTEN)
#include <emscripten.h>
#include <emscripten/bind.h>
#endif

#include "graphics.hpp"
#include "shader.hpp"
#include "entities.hpp"
#include "loader.hpp"
#include "controls.hpp"

GLFWwindow *window;
Entities entities;
Camera camera;

#if defined(EMSCRIPTEN)

void emLoadPdbFile(std::string path) {
    entities.clear();

    // Load default file
    PdbFile pdb_file;
    loadPdbFile(pdb_file, path);
    createEntitiesFromPdbFile(entities, pdb_file, camera);
}
void emLoadMolFile(std::string path) {
    entities.clear();

    // Load default file
    MolFile mol_file;
    loadMolFile(mol_file, path);
    createEntitiesFromMolFile(entities, mol_file, camera);
}
// Binding c++ functions to js api
EMSCRIPTEN_BINDINGS(module)
{
    emscripten::function("loadPdbFile", &emLoadPdbFile);
    emscripten::function("loadMolFile", &emLoadMolFile);
}
#endif

void terminate()
{
    glfwTerminate();
}

bool initialize()
{
    // Initialize the glfw library
    //glfwSetErrorCallback(error_callback);
    if(!glfwInit()) return false;
        
    // Output glfw version.
    std::cout << glfwGetVersionString() << "\n";

    // Create main window.
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, 1);
    glfwWindowHint(GLFW_MAXIMIZED, 1);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, 1);
#if defined(EMSCRIPTEN)
    shader::glsl_version = "#version 300 es\n";
#else 
    shader::glsl_version = "#version 330 core\n";
#endif

#if defined(EMSCRIPTEN)
    double w,h;
    emscripten_get_element_css_size("#canvas", &w, &h);
    window_width = (int)w;
    window_height = (int)h; 
    emscripten_set_canvas_element_size("#canvas", window_width, window_height);
#else
    window_width = 640;
    window_height = 480;
#endif 

    window = glfwCreateWindow(window_width, window_height, "app [glfw]", 0, 0);
    if(window == NULL) {
        std::cerr << "Failed to create glfw window\n";
        terminate();
        return false;
    }
    glfwMakeContextCurrent(window); 

#if defined(EMSCRIPTEN)
    // 
    // Setup callbacks
    //
    EMSCRIPTEN_RESULT res = emscripten_set_resize_callback(EMSCRIPTEN_EVENT_TARGET_WINDOW, 0 /*void *userData*/, true /* useCapture */, emResizeCallback);
    res |= emscripten_set_wheel_callback( "#canvas", 0 /*void *userData*/, true /* useCapture */, emScrollCallback);
    res |= emscripten_set_mousedown_callback( "#canvas", 0 /*void *userData*/, true /* useCapture */, emMouseDownCallback);
    res |= emscripten_set_mousemove_callback( "#canvas", 0 /*void *userData*/, true /* useCapture */, emMouseMoveCallback);
    res |= emscripten_set_mouseup_callback( "#canvas", 0 /*void *userData*/, true /* useCapture */, emMouseUpCallback);

    res |= emscripten_set_touchstart_callback( "#canvas", 0 /*void *userData*/, true /* useCapture */, emTouchStartCallback);
    res |= emscripten_set_touchmove_callback( "#canvas", 0 /*void *userData*/, true /* useCapture */, emTouchMoveCallback);
    res |= emscripten_set_touchend_callback( "#canvas", 0 /*void *userData*/, true /* useCapture */, emTouchEndCallback);
    res |= emscripten_set_touchcancel_callback( "#canvas", 0 /*void *userData*/, true /* useCapture */, emTouchCancelCallback);

#else
    auto primary_monitor = glfwGetPrimaryMonitor();
    GLFWmonitor *monitor = glfwGetWindowMonitor(window);
    GLFWvidmode const *mode = glfwGetVideoMode(monitor ? monitor : primary_monitor);
    if (monitor)
        glfwSetWindowMonitor(window, 0, 0, 0, window_width, window_height, mode->refreshRate);
    else
    {
        glfwSetWindowMonitor(window, primary_monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
    }

    glfwSetWindowSizeCallback(window, &glfwResizeCallback);
    glfwSetScrollCallback(window, &glfwScrollCallback);
#endif

    //
    // Init OpenGL
    //
#if !defined(EMSCRIPTEN)
    auto ret = gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
#endif

#if defined(EMSCRIPTEN)
    loadBasicShader("data/shaders/basic.gl");
    loadBasicInstancedShader("data/shaders/basic_instanced.gl");
#else
    // @todo Live update
    loadBasicShader("data/shaders/basic.gl");
    loadBasicInstancedShader("data/shaders/basic_instanced.gl");
#endif

    createDefaultCamera(camera);
    initControls(window);
    initGraphics();

    // Load default file
    PdbFile pdb_file;
    loadPdbFile(pdb_file, "data/examples/pdb/1bzv.pdb");
    createEntitiesFromPdbFile(entities, pdb_file, camera);

    return (glGetError() == GL_NO_ERROR);
}

void step()
{
    static double prev_time = glfwGetTime();
    double curr_time = glfwGetTime();

    double delta = curr_time - prev_time;

#if defined(EMSCRIPTEN)
    if(delta <= 1.0/60.0) 
        return;
#endif

    handleControls(window, camera, delta);
    if(window_resized) {
        updateCamera(camera);
        glfwSetWindowSize(window, window_width, window_height);
    }

    // Rendering
    {
        glViewport(0, 0, window_width, window_height);

        glClearColor(0.5, 1.0, 0.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT);
        bindBackbuffer();
        drawEntities(entities, camera);
    }
    // swap backbuffer with front buffer
    glfwSwapBuffers(window);
    window_resized = false;
    glfwPollEvents();

#if DEBUG
    GLenum code = glGetError();
    if(code != GL_NO_ERROR){
        static const GLubyte* err;
        err = gluErrorString(code);
        std::cerr <<  "<--------------------OpenGL ERROR-------------------->\n" << err << "\n";
    }
#endif

    prev_time = curr_time;
}

void runMessageLoop()
{
#if defined(EMSCRIPTEN)
    emscripten_set_main_loop(step, 0 /* fps */, 1 /* simulate_infinite_loop */);
#else
    while (!glfwWindowShouldClose(window))
    {
        step();
    }
    terminate();
#endif
}

int main() {
    if (initialize()) {
        runMessageLoop();
    } else {
        return 1;
    }

    return 0;
}

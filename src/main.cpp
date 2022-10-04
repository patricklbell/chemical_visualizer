#include <array>
#include <chrono>
#include <iostream>
#include <map>
#include <thread>

#if defined(EMSCRIPTEN)
#include <emscripten.h>
#include <emscripten/bind.h>
#endif

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include "graphics.hpp"
#include "shader.hpp"
#include "entities.hpp"
#include "loader.hpp"
#include "controls.hpp"
#include "ui.hpp"

Camera camera;
Entities entities;

#if defined(EMSCRIPTEN)
void emLoadPdbFileHelper(std::string path) {
    loadPdbFileHelper(path, entities, camera);
}
void emLoadMolFileHelper(std::string path) {
    loadMolFileHelper(path, entities, camera);
}
UiMode emGetUiMode() {
    return ui::mode;
}
PdbDrawSettings emGetPdbDrawSettings() {
    return ui::pdbfile_settings;
}
void emSetPdbDrawSettings(PdbDrawSettings settings) {
    if(ui::mode == UiMode::PDB && ui::pdbfile.models.size() > 0 && ui::pdbfile.models[0].renderable) {
        bool entities_changed = (settings.draw_hetero_atoms != ui::pdbfile_settings.draw_hetero_atoms) || 
                                (settings.draw_water_atoms != ui::pdbfile_settings.draw_water_atoms) || 
                                (settings.draw_residue_atoms != ui::pdbfile_settings.draw_residue_atoms) || 
                                (settings.draw_residue_ribbons != ui::pdbfile_settings.draw_residue_ribbons);
        bool colors_changed   = (settings.hetero_atoms_alpha != ui::pdbfile_settings.hetero_atoms_alpha) || 
                                (settings.water_atoms_alpha != ui::pdbfile_settings.water_atoms_alpha) || 
                                (settings.residue_atoms_alpha != ui::pdbfile_settings.residue_atoms_alpha) || 
                                (settings.residue_ribbons_alpha != ui::pdbfile_settings.residue_ribbons_alpha) || 
                                (settings.residue_color_mode != ui::pdbfile_settings.residue_color_mode);

        if (entities_changed) {
            entities.clear();
            createEntitiesFromPdbModel(entities, ui::pdbfile.models[ui::selected_pdb_model_serial], settings, camera);
        } else if (colors_changed) {
            updateEntityColorsFromPdbModel(entities, ui::pdbfile.models[ui::selected_pdb_model_serial], settings);
        }
    }

    ui::pdbfile_settings = settings;
} 
void emSetPerspectiveCamera() {
    camera.is_ortho = false;
    updateCamera(camera);
}
void emSetOrthographicCamera() {
    camera.is_ortho = true;
    updateCameraProjection(camera);
}
void emSetLightMode() {
    ui::dark_mode = false;
}
void emSetDarkMode() {
    ui::dark_mode = true;
}
void emSetFullscreen() {
    ui::fullscreen = true;
}
void emUnsetFullscreen() {
    ui::fullscreen = false;
}

// Binding c++ functions to js api
EMSCRIPTEN_BINDINGS(module)
{
    // Define datatype bindings
    emscripten::enum_<UiMode>("UiMode")
        .value("NONE", UiMode::NONE)
        .value("PDB", UiMode::PDB)
        .value("MOL", UiMode::MOL)
        ;
    emscripten::enum_<PdbResidueColorMode>("PdbResidueColorMode")
        .value("CHAIN", PdbResidueColorMode::CHAIN)
        .value("SECONDARY", PdbResidueColorMode::SECONDARY)
        .value("AMINO_ACID", PdbResidueColorMode::AMINO_ACID)
        ;
    emscripten::value_object<PdbDrawSettings>("PdbDrawSettings")
        .field("draw_hetero_atoms", &PdbDrawSettings::draw_hetero_atoms)
        .field("draw_water_atoms", &PdbDrawSettings::draw_water_atoms)
        .field("draw_residue_atoms", &PdbDrawSettings::draw_residue_atoms)
        .field("draw_residue_ribbons", &PdbDrawSettings::draw_residue_ribbons)
        .field("hetero_atoms_alpha", &PdbDrawSettings::hetero_atoms_alpha)
        .field("water_atoms_alpha", &PdbDrawSettings::water_atoms_alpha)
        .field("residue_atoms_alpha", &PdbDrawSettings::residue_atoms_alpha)
        .field("residue_ribbons_alpha", &PdbDrawSettings::residue_ribbons_alpha)
        .field("residue_color_mode", &PdbDrawSettings::residue_color_mode)
        ;

    emscripten::function("loadPdbFile", &emLoadPdbFileHelper);
    emscripten::function("loadMolFile", &emLoadMolFileHelper);
    emscripten::function("getUiMode", &emGetUiMode);
    emscripten::function("getPdbDrawSettings", &emGetPdbDrawSettings);
    emscripten::function("setPdbDrawSettings", &emSetPdbDrawSettings);
    emscripten::function("setPerspectiveCamera", &emSetPerspectiveCamera);
    emscripten::function("setOrthographicCamera", &emSetOrthographicCamera);
    emscripten::function("setLightMode", &emSetLightMode);
    emscripten::function("setDarkMode", &emSetDarkMode);
    emscripten::function("setFullscreen", &emSetFullscreen);
    emscripten::function("unsetFullscreen", &emUnsetFullscreen);
}
#endif

void cleanup()
{
    deleteShaderPrograms();
    glfwTerminate();
}

bool initialize()
{
    // Initialize the glfw library
    //glfwSetErrorCallback(error_callback);
    if(!glfwInit()) return false;

    // Create main window.
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, 1);
#if defined(EMSCRIPTEN)
    shader::glsl_version = "#version 300 es\n#define GLES\n";
#else 
    shader::glsl_version = "#version 330 core\n";
#endif

#if defined(EMSCRIPTEN)
    glfwWindowHint(GLFW_MAXIMIZED, 1);

    double w,h;
    emscripten_get_element_css_size("#canvas", &w, &h);
    window_width = (int)w;
    window_height = (int)h; 
#else
    glfwWindowHint(GLFW_RESIZABLE, 1);
    window_width = 1024;
    window_height = 700;
#endif 

    window = glfwCreateWindow(window_width, window_height, "Chemical Visualizer", NULL, NULL);
    if(window == NULL) {
        fprintf(stderr, "Failed to create glfw window\n");
        cleanup();
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
    glfwSetWindowSizeCallback(window, &glfwResizeCallback);
    glfwSetScrollCallback(window, &glfwScrollCallback);
#endif

    //
    // Init OpenGL
    //
#if !defined(EMSCRIPTEN)
    if (glewInit() != GLEW_OK)
    {
        fprintf(stderr, "Failed to initialize GLEW\n");
        glfwTerminate();
        return false;
    }
    glEnable(GL_MULTISAMPLE);
#endif

    loadBasicShader("data/shaders/basic.gl");
    loadBasicInstancedShader("data/shaders/basic_instanced.gl");

    createDefaultCamera(camera);
    initControls(window);
    initGraphics();
#if !defined(EMSCRIPTEN)
    initGui();
#endif

    // Cache this so we don't have to load whole file, then it can be used with emscripten
#if !defined(EMSCRIPTEN)
    // Not necessary if you don't want to see peptide's actual bonds
    loadPdbDictionaryFile(ui::pdb_dictionary, "data/examples/pdb/het_dictionary.pdb");
#endif

    // Load default file
    {
        loadPdbFileHelper("data/examples/pdb/1bzv.pdb", entities, camera);
    };

    return true;
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
        window_resized = false;
    }

    // Rendering
    {
        glViewport(0, 0, window_width, window_height);

        glClearColor(0.5, 1.0, 0.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT);
        bindBackbuffer();
        drawEntities(entities, camera);

#if !defined(EMSCRIPTEN)
        checkWriteFrambufferToTga();
        drawGui(camera, entities);
#endif
    }
    // swap backbuffer with front buffer
    glfwSwapBuffers(window);
    glfwPollEvents();

#if DEBUG
    GLenum code = glGetError();
    if(code != GL_NO_ERROR){
        static const GLubyte* err;
        err = gluErrorString(code);
        fprintf(stderr, "<--------------------OpenGL ERROR-------------------->\n%s\n", err);
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
    cleanup();
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

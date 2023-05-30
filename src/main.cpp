#ifdef EMSCRIPTEN
#include <emscripten.h>
#include <emscripten/bind.h>
#endif

#include <string>
#include <iostream>

#include <platform.hpp>
#include <globals.hpp>
#include <graphics.hpp>
#include <controls.hpp>
#include <utilities/strings.hpp>
#include <shader/globals.hpp>

namespace Globals {
    GLFWwindow* window;

    Camera camera              = Camera();
    Entities entities          = Entities();
    DrawSettings draw_settings = DrawSettings();

    MolFile mol;

    PdbDictionary pdb_dictionary;
    PdbDrawSettings pdb_settings;
    int pdb_selected_model = 0;
    PdbFile pdb;
}    // namespace Globals
using namespace Globals;

void loop() {
    static double ptime = glfwGetTime();
    double ctime        = glfwGetTime();
    double dt           = ctime - ptime;

#ifdef EMSCRIPTEN
    if (dt <= 1.0 / 60.0)
        return;
#endif

    handle_controls(camera, dt);
    camera.update();

    draw_entities(entities, camera, draw_settings);

    glfwSwapBuffers(window);
    glfwPollEvents();

#if !defined(EMSCRIPTEN) && !defined(NDEBUG)
    live_update_shaders();

    GLenum code = glGetError();
    if (code != GL_NO_ERROR)
        std::cerr << "OpenGL Error: " << gluErrorString(code) << "\n";
#endif

    ptime = ctime;
}

int main(int argc, char* argv[]) {
#ifndef EMSCRIPTEN
    const char* msg =
        "cviz: cviz [FILE]\n"
        "    Visualize the 3D structure of a .pdb or .mol FILE\n";
    if (argc < 2) {
        std::cout << msg;
        return 1;
    }
#endif

    if (!init())
        return 1;

#if EMSCRIPTEN
    emscripten_set_main_loop(loop, 0, true);
#else
    std::string path(argv[1]);
    if (endsWith(path, ".pdb")) {
        loadPdbFile(pdb, path, pdb_dictionary);

        if (pdb.models.size() > 0) {
            createPdbModelMeshes(pdb.models[pdb_selected_model]);
            createEntitiesFromPdbModel(entities, pdb.models[pdb_selected_model], pdb_settings, camera);
        }
    } else if (endsWith(path, ".mol")) {
        loadMolFile(mol, path);
        createEntitiesFromMolFile(entities, mol, camera);
    } else {
        std::cout << msg;
        return 1;
    }

    while (!glfwWindowShouldClose(window))
        loop();
    glfwTerminate();
#endif

    return 0;
}

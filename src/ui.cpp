#include <filesystem>
#include <stack>
#include <limits>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>

// Include ImGui
#include "glm/detail/func_geometric.hpp"
#include "glm/detail/type_mat.hpp"
#include "glm/fwd.hpp"
#include "glm/gtc/quaternion.hpp"
#include "imgui.h"
#include "ImGuizmo.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imfilebrowser.hpp"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <string>

#include "ui.hpp"
#include "globals.hpp"
#include "utilities.hpp"
#include "graphics.hpp"
#include "shader.hpp"
#include "assets.hpp"
#include "controls.hpp"
#include "entities.hpp"
#include "loader.hpp"

namespace editor {
    std::string im_file_dialog_type;
    ImGui::FileBrowser im_file_dialog(ImGuiFileBrowserFlags_EnterNewFilename | ImGuiFileBrowserFlags_NoTitleBar);
}
using namespace editor;

void initGui(){
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    auto &io = ImGui::GetIO();
    (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    
    // create a file browser instance
    im_file_dialog_type = "";

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(shader::glsl_version.c_str());
}

void drawGui(Camera &camera, Entities &entities){
    // Start the Dear ImGui frame;
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    
    {
        ImGui::SetNextWindowPos(ImVec2(0,0));
        ImGui::SetNextWindowSizeConstraints(ImVec2(200, window_height), ImVec2(window_width / 2.0, window_height));
        ImGui::Begin("Load Files", NULL, ImGuiWindowFlags_NoMove);

        if (ImGui::Button("Load Molecule File", ImVec2(ImGui::GetWindowWidth()-10, 20))){
            im_file_dialog_type = "loadMol";
            im_file_dialog.SetCurrentTypeFilterIndex(2);
            im_file_dialog.SetTypeFilters({".mol"});

            im_file_dialog.SetPwd("data/examples/molfiles");
            im_file_dialog.Open();
        }
        if (ImGui::Button("Load Protein Database File", ImVec2(ImGui::GetWindowWidth()-10, 20))){
            im_file_dialog_type = "loadPdb";
            im_file_dialog.SetCurrentTypeFilterIndex(1);
            im_file_dialog.SetTypeFilters({".pdb"});

            im_file_dialog.SetPwd("data/examples/pdb");
            im_file_dialog.Open();
        }


        if (ImGui::Button("Print Camera", ImVec2(ImGui::GetWindowWidth()-10, 20))){
            printf("Position: %f, %f, %f, Target: %f, %f, %f\n", camera.position.x, camera.position.y, camera.position.z, camera.target.x, camera.target.y, camera.target.z);
        }

        if (ImGui::Button("Save Frame", ImVec2(ImGui::GetWindowWidth()-10, 20))){
            auto t = std::time(nullptr);
            auto tm = *std::localtime(&t);

            std::ostringstream oss;
            oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
            auto path = "data/screenshots/" + oss.str() + ".tga";

            pushWriteFramebufferToTga(path);
        }

        ImGui::SetCursorPosY(window_height - ImGui::GetTextLineHeightWithSpacing()*3);
        ImGui::TextWrapped("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    }
    ImGui::End();

    im_file_dialog.Display();
    if(im_file_dialog.HasSelected()){
        auto p = im_file_dialog.GetSelected().generic_string();
        printf("Selected filename: %s\n", im_file_dialog.GetSelected().c_str());
        if(im_file_dialog_type == "loadMol"){
            MolFile molfile;
            loadMolFile(molfile, p);

            entities.clear();
            createEntitiesFromMolFile(entities, molfile, camera);
        } else if(im_file_dialog_type == "loadPdb"){
            PdbFile pdbfile;
            loadPdbFile(pdbfile, p, &pdb_dictionary);

            entities.clear();
            createEntitiesFromPdbFile(entities, pdbfile, camera);
        } else {
            fprintf(stderr, "Unhandled imgui file dialog type %s.\n", p.c_str());
        }
        im_file_dialog.ClearSelected();
    }

    ImGui::Render();
    auto &io = ImGui::GetIO();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

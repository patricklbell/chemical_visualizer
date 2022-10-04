#include <filesystem>
#include <stack>
#include <limits>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

// Include ImGui
#include "glm/detail/func_geometric.hpp"
#include "glm/detail/type_mat.hpp"
#include "glm/fwd.hpp"
#include "glm/gtc/quaternion.hpp"

#if !defined(EMSCRIPTEN)
#include "imgui.h"
#include "ImGuizmo.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imfilebrowser.hpp"
#endif

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <string>

#include "ui.hpp"
#include "utilities.hpp"
#include "graphics.hpp"
#include "shader.hpp"
#include "assets.hpp"
#include "controls.hpp"
#include "entities.hpp"
#include "loader.hpp"

namespace ui {
#if defined(EMSCRIPTEN)
    bool fullscreen = false;
#else
    std::string im_file_dialog_type;
    ImGui::FileBrowser im_file_dialog(ImGuiFileBrowserFlags_EnterNewFilename | ImGuiFileBrowserFlags_NoTitleBar);
#endif

    UiMode mode = UiMode::NONE;

    MolFile molfile;

    PdbFile pdbfile;
    PdbDictionary pdb_dictionary;
    PdbDrawSettings pdbfile_settings;
    int selected_pdb_model_serial = 0;

    std::string loaded_file_path = "";

    bool dark_mode = true;
}
using namespace ui;

void loadPdbFileHelper(std::string path, Entities &entities, Camera &camera) {
    ui::pdbfile.clear();
    loadPdbFile(ui::pdbfile, path, &pdb_dictionary);

    if (ui::pdbfile.models.size() > 0) {
        entities.clear();
        createPdbModelMeshes(ui::pdbfile.models[0]);
        createEntitiesFromPdbModel(entities, ui::pdbfile.models[0], ui::pdbfile_settings, camera);

        positionCameraSelectionInView(camera, 0.7);

        ui::loaded_file_path = path;
        ui::mode = UiMode::PDB;
    }
}
void loadMolFileHelper(std::string path, Entities &entities, Camera &camera) {
    molfile.clear();
    loadMolFile(molfile, path);

    entities.clear();
    createEntitiesFromMolFile(entities, molfile, camera);

    positionCameraSelectionInView(camera, 0.8);

    loaded_file_path = path;
    mode = UiMode::MOL;
}

#if !defined(EMSCRIPTEN)
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

    auto& style = ImGui::GetStyle();

    // Setup Dear ImGui style
    style.FrameRounding = 4.0f;
    style.GrabRounding = 4.0f;
    style.FramePadding = ImVec2(4.0, 4.0);
    style.ItemSpacing = ImVec2(4.0, 8.0);
    style.ItemInnerSpacing = ImVec2(5.0, 1.0);
    style.CellPadding = ImVec2(1.0, 6.0);
    style.Alpha = 0.95;

    ImVec4* colors = style.Colors;
    colors[ImGuiCol_Text] = ImVec4(0.95f, 0.96f, 0.98f, 1.00f);
    colors[ImGuiCol_TextDisabled] = ImVec4(0.36f, 0.42f, 0.47f, 1.00f);
    colors[ImGuiCol_WindowBg] = ImVec4(0.11f, 0.15f, 0.17f, 1.00f);
    colors[ImGuiCol_ChildBg] = ImVec4(0.15f, 0.18f, 0.22f, 1.00f);
    colors[ImGuiCol_PopupBg] = ImVec4(0.08f, 0.08f, 0.08f, 0.94f);
    colors[ImGuiCol_Border] = ImVec4(0.08f, 0.10f, 0.12f, 1.00f);
    colors[ImGuiCol_BorderShadow] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_FrameBg] = ImVec4(0.20f, 0.25f, 0.29f, 1.00f);
    colors[ImGuiCol_FrameBgHovered] = ImVec4(0.12f, 0.20f, 0.28f, 1.00f);
    colors[ImGuiCol_FrameBgActive] = ImVec4(0.09f, 0.12f, 0.14f, 1.00f);
    colors[ImGuiCol_TitleBg] = ImVec4(0.09f, 0.12f, 0.14f, 0.65f);
    colors[ImGuiCol_TitleBgActive] = ImVec4(0.08f, 0.10f, 0.12f, 1.00f);
    colors[ImGuiCol_TitleBgCollapsed] = ImVec4(0.00f, 0.00f, 0.00f, 0.51f);
    colors[ImGuiCol_MenuBarBg] = ImVec4(0.15f, 0.18f, 0.22f, 1.00f);
    colors[ImGuiCol_ScrollbarBg] = ImVec4(0.02f, 0.02f, 0.02f, 0.39f);
    colors[ImGuiCol_ScrollbarGrab] = ImVec4(0.20f, 0.25f, 0.29f, 1.00f);
    colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.18f, 0.22f, 0.25f, 1.00f);
    colors[ImGuiCol_ScrollbarGrabActive] = ImVec4(0.09f, 0.21f, 0.31f, 1.00f);
    colors[ImGuiCol_CheckMark] = ImVec4(0.28f, 0.56f, 1.00f, 1.00f);
    colors[ImGuiCol_SliderGrab] = ImVec4(0.28f, 0.56f, 1.00f, 1.00f);
    colors[ImGuiCol_SliderGrabActive] = ImVec4(0.37f, 0.61f, 1.00f, 1.00f);
    colors[ImGuiCol_Button] = ImVec4(0.20f, 0.25f, 0.29f, 1.00f);
    colors[ImGuiCol_ButtonHovered] = ImVec4(0.28f, 0.56f, 1.00f, 1.00f);
    colors[ImGuiCol_ButtonActive] = ImVec4(0.06f, 0.53f, 0.98f, 1.00f);
    colors[ImGuiCol_Header] = ImVec4(0.20f, 0.25f, 0.29f, 0.55f);
    colors[ImGuiCol_HeaderHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.80f);
    colors[ImGuiCol_HeaderActive] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
    colors[ImGuiCol_Separator] = ImVec4(0.20f, 0.25f, 0.29f, 1.00f);
    colors[ImGuiCol_SeparatorHovered] = ImVec4(0.10f, 0.40f, 0.75f, 0.78f);
    colors[ImGuiCol_SeparatorActive] = ImVec4(0.10f, 0.40f, 0.75f, 1.00f);
    colors[ImGuiCol_ResizeGrip] = ImVec4(0.26f, 0.59f, 0.98f, 0.25f);
    colors[ImGuiCol_ResizeGripHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.67f);
    colors[ImGuiCol_ResizeGripActive] = ImVec4(0.26f, 0.59f, 0.98f, 0.95f);
    colors[ImGuiCol_Tab] = ImVec4(0.11f, 0.15f, 0.17f, 1.00f);
    colors[ImGuiCol_TabHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.80f);
    colors[ImGuiCol_TabActive] = ImVec4(0.20f, 0.25f, 0.29f, 1.00f);
    colors[ImGuiCol_TabUnfocused] = ImVec4(0.11f, 0.15f, 0.17f, 1.00f);
    colors[ImGuiCol_TabUnfocusedActive] = ImVec4(0.11f, 0.15f, 0.17f, 1.00f);
    colors[ImGuiCol_PlotLines] = ImVec4(0.61f, 0.61f, 0.61f, 1.00f);
    colors[ImGuiCol_PlotLinesHovered] = ImVec4(1.00f, 0.43f, 0.35f, 1.00f);
    colors[ImGuiCol_PlotHistogram] = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
    colors[ImGuiCol_PlotHistogramHovered] = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
    colors[ImGuiCol_TextSelectedBg] = ImVec4(0.26f, 0.59f, 0.98f, 0.35f);
    colors[ImGuiCol_DragDropTarget] = ImVec4(1.00f, 1.00f, 0.00f, 0.90f);
    colors[ImGuiCol_NavHighlight] = ImVec4(0.26f, 0.59f, 0.98f, 1.00f);
    colors[ImGuiCol_NavWindowingHighlight] = ImVec4(1.00f, 1.00f, 1.00f, 0.70f);
    colors[ImGuiCol_NavWindowingDimBg] = ImVec4(0.80f, 0.80f, 0.80f, 0.20f);
    colors[ImGuiCol_ModalWindowDimBg] = ImVec4(0.80f, 0.80f, 0.80f, 0.35f);

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
        ImGui::SetNextWindowSizeConstraints(ImVec2(170, window_height), ImVec2(window_width / 2.0, window_height));
        ImGui::SetNextWindowSize(ImVec2(250, window_height), ImGuiCond_FirstUseEver);
        ImGui::Begin("###main-panel", NULL, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoTitleBar);

        auto button_size = ImVec2(ImGui::GetWindowWidth() - 10, 20);
        if (ImGui::Button("Load Molecule File", button_size)){
            im_file_dialog_type = "loadMol";
            im_file_dialog.SetCurrentTypeFilterIndex(2);
            im_file_dialog.SetTypeFilters({".mol"});

            im_file_dialog.SetPwd("data/examples/molfiles");
            im_file_dialog.Open();
        }
        if (ImGui::Button("Load Protein Database File", button_size)){
            im_file_dialog_type = "loadPdb";
            im_file_dialog.SetCurrentTypeFilterIndex(1);
            im_file_dialog.SetTypeFilters({".pdb"});

            im_file_dialog.SetPwd("data/examples/pdb");
            im_file_dialog.Open();
        }

        ImGui::SetCursorPosY(ImGui::GetCursorPosY() + ImGui::GetTextLineHeightWithSpacing() * 0.5);
        if (ImGui::CollapsingHeader("Camera")) {
            if (ImGui::Button("Save Frame", button_size)) {
                auto t = std::time(nullptr);
                auto tm = *std::localtime(&t);

                std::ostringstream oss;
                oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
                auto path = "data/screenshots/" + oss.str() + ".tga";

                pushWriteFramebufferToTga(path);
            }
            if(ImGui::Checkbox("Orthographic", &camera.is_ortho)) {
                updateCameraProjection(camera);
            }
            ImGui::Checkbox("Dark Background", &dark_mode);
        }

        switch (mode)
        {
        case UiMode::NONE:
            break;
        case UiMode::PDB:
            ImGui::SetCursorPosY(ImGui::GetCursorPosY() + ImGui::GetTextLineHeightWithSpacing() * 0.5);
            if (ImGui::CollapsingHeader("Pdb Options", ImGuiTreeNodeFlags_DefaultOpen)) {
                bool colors_changed = false;
                bool entities_changed = false;

                // @todo load model names from pdb
                if (pdbfile.models.size() > 1) {
                    std::string model_name = std::to_string(pdbfile.models[selected_pdb_model_serial].serial);

                    ImGui::SetCursorPosY(ImGui::GetCursorPosY() + ImGui::GetTextLineHeightWithSpacing() * 0.5);
                    ImGui::Text("List of Models");
                    ImGui::SetNextItemWidth(button_size.x);
                    if (ImGui::BeginCombo("##model-combo", model_name.c_str())) {
                        for(auto &p : pdbfile.models) {
                            auto& m = p.second;

                            bool is_selected = (m.serial == selected_pdb_model_serial);
                            if (is_selected)
                                ImGui::SetItemDefaultFocus();

                            if (ImGui::Selectable(std::to_string(m.serial).c_str(), is_selected)) {
                                selected_pdb_model_serial = m.serial;
                                if (!m.renderable) {
                                    createPdbModelMeshes(m);
                                }

                                entities.clear();
                                createEntitiesFromPdbModel(entities, m, pdbfile_settings, camera);
                            }
                        }
                        ImGui::EndCombo();
                    }
                }
                entities_changed |= ImGui::Checkbox("Heterogen Molecules", &pdbfile_settings.draw_hetero_atoms);
                entities_changed |= ImGui::Checkbox("Water Molecules", &pdbfile_settings.draw_water_atoms);
                entities_changed |= ImGui::Checkbox("Residue Molecules", &pdbfile_settings.draw_residue_atoms);
                entities_changed |= ImGui::Checkbox("Residue Ribbon", &pdbfile_settings.draw_residue_ribbons);

                ImGui::SetCursorPosY(ImGui::GetCursorPosY() + ImGui::GetTextLineHeightWithSpacing() * 0.5);
                if (ImGui::CollapsingHeader("Color Options")) {
                    const std::string color_modes[PdbResidueColorMode::NUM_MODES] = { "Chain", "Secondary", "Amino Acids" };
                    std::string mode_name = color_modes[(int)pdbfile_settings.residue_color_mode];

                    ImGui::Text("Residue Coloring");
                    ImGui::SetNextItemWidth(button_size.x);
                    if (ImGui::BeginCombo("##residue-color-combo", mode_name.c_str())) {
                        for(int i = 0; i < (int)PdbResidueColorMode::NUM_MODES; ++i) {
                            bool is_selected = ((int)pdbfile_settings.residue_color_mode == i);
                            if (is_selected)
                                ImGui::SetItemDefaultFocus();

                            if (ImGui::Selectable(color_modes[i].c_str(), is_selected)) {
                                pdbfile_settings.residue_color_mode = (PdbResidueColorMode)i;
                                colors_changed = true;
                            }
                        }
                        ImGui::EndCombo();
                    }

                    ImGui::Text("Transparency");
                    if (pdbfile_settings.draw_hetero_atoms   ) colors_changed |= ImGui::SliderFloat("Heterogen",    &pdbfile_settings.hetero_atoms_alpha   , 0.0f, 1.0f, "%.2f");
                    if (pdbfile_settings.draw_water_atoms    ) colors_changed |= ImGui::SliderFloat("Water",        &pdbfile_settings.water_atoms_alpha    , 0.0f, 1.0f, "%.2f");
                    if (pdbfile_settings.draw_residue_atoms  ) colors_changed |= ImGui::SliderFloat("Residue",      &pdbfile_settings.residue_atoms_alpha  , 0.0f, 1.0f, "%.2f");
                    if (pdbfile_settings.draw_residue_ribbons) colors_changed |= ImGui::SliderFloat("Ribbon",       &pdbfile_settings.residue_ribbons_alpha, 0.0f, 1.0f, "%.2f");
                }

                if (entities_changed) {
                    entities.clear();
                    createEntitiesFromPdbModel(entities, pdbfile.models[selected_pdb_model_serial], pdbfile_settings, camera);
                } else if (colors_changed) {
                    updateEntityColorsFromPdbModel(entities, pdbfile.models[selected_pdb_model_serial], pdbfile_settings);
                }
            }
            break;
        case UiMode::MOL:
            break;
        default:
            break;
        }

        ImGui::SetCursorPosY(window_height - ImGui::GetTextLineHeightWithSpacing()*2);
        ImGui::Text("%.3f ms/f, %.1f FPS", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
        std::string displayed_file_path;
        if (loaded_file_path.size() > 16) {
            displayed_file_path = loaded_file_path.substr(0, 4) + ".." + loaded_file_path.substr(loaded_file_path.size() - 12, 12);
        }
        else {
            displayed_file_path = loaded_file_path;
        }
        if(mode != UiMode::NONE) ImGui::Text("Viewing file %s", displayed_file_path.c_str());
    }
    ImGui::End();

    im_file_dialog.Display();
    if(im_file_dialog.HasSelected()){
        auto p = im_file_dialog.GetSelected().generic_string();
        printf("Selected filename: %s\n", im_file_dialog.GetSelected().c_str());
        if(im_file_dialog_type == "loadMol") {
            loadMolFileHelper(p, entities, camera);
        } else if(im_file_dialog_type == "loadPdb"){
            loadPdbFileHelper(p, entities, camera);
        } else {
            fprintf(stderr, "Unhandled imgui file dialog type %s.\n", p.c_str());
        }
        im_file_dialog.ClearSelected();
    }

    ImGui::Render();
    auto &io = ImGui::GetIO();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
#endif
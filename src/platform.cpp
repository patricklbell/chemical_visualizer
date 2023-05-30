#include <platform.hpp>

#include <array>
#include <chrono>
#include <iostream>
#include <sstream>
#include <map>
#include <thread>

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include <graphics.hpp>
#include <globals.hpp>
#include <controls.hpp>
#include <loader.hpp>
#include <shader/globals.hpp>

using namespace Globals;

#ifdef EMSCRIPTEN
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include <emscripten/html5.h>
PdbDrawSettings emGetPdbDrawSettings() { return pdb_settings; }
void emSetPdbDrawSettings(PdbDrawSettings settings) {
    bool entities_changed = (settings.draw_hetero_atoms != pdb_settings.draw_hetero_atoms) ||
                            (settings.draw_water_atoms != pdb_settings.draw_water_atoms) ||
                            (settings.draw_residue_atoms != pdb_settings.draw_residue_atoms) ||
                            (settings.draw_residue_ribbons != pdb_settings.draw_residue_ribbons);
    bool colors_changed = (settings.hetero_atoms_alpha != pdb_settings.hetero_atoms_alpha) ||
                          (settings.water_atoms_alpha != pdb_settings.water_atoms_alpha) ||
                          (settings.residue_atoms_alpha != pdb_settings.residue_atoms_alpha) ||
                          (settings.residue_ribbons_alpha != pdb_settings.residue_ribbons_alpha) ||
                          (settings.residue_color_mode != pdb_settings.residue_color_mode);

    pdb_settings = settings;

    if (!pdb.models.size() || !pdb.models[0].renderable)
        return;

    if (entities_changed) {
        entities.clear();
        auto camera_state = camera;
        createEntitiesFromPdbModel(entities, pdb.models[pdb_selected_model], settings,
                                   camera);
        camera = camera_state;
    } else if (colors_changed) {
        updateEntityColorsFromPdbModel(entities, pdb.models[pdb_selected_model], settings);
    }
}

Camera::Type emGetCameraType() { return camera.type; }
void emSetCameraType(Camera::Type type) { camera.set_type(type); }

DrawSettings emGetDrawSettings() { return draw_settings; }
void emSetDrawSettings(DrawSettings _draw_settings) { draw_settings = _draw_settings; }

void emSetViewportSize(int w, int h) {
    window_width  = w;
    window_height = h;

    camera.set_aspect_ratio((double)window_width / (double)window_height);
    glfwSetWindowSize(window, window_width, window_height);
}

void emLoadPdb(std::string contents) {
    pdb.clear();
    mol.clear();
    entities.clear();

    auto stream = std::istringstream(contents);
    loadPdbStream(pdb, stream, pdb_dictionary);
    if (pdb.models.size()) {
        pdb_selected_model = 0;
        createPdbModelMeshes(pdb.models[pdb_selected_model]);
        createEntitiesFromPdbModel(entities, pdb.models[pdb_selected_model], pdb_settings, camera);
    }
}

void emLoadMol(std::string contents) {
    pdb.clear();
    mol.clear();
    entities.clear();

    auto stream = std::istringstream(contents);
    loadMolStream(mol, stream);
    createEntitiesFromMolFile(entities, mol, camera);
}

void emLoadPdbDictionary(std::string contents) {
    auto stream = std::istringstream(contents);
    loadPdbDictionaryStream(pdb_dictionary, stream);

    if (pdb.models.size()) {
        entities.clear();
        addDictionaryToPdb(pdb, pdb_dictionary);
        createEntitiesFromPdbModel(entities, pdb.models[pdb_selected_model], pdb_settings, camera);
    }
}

struct EmPdbResidue {
    std::string id;
    std::string name;
};
struct EmPdbChain {
    std::string id;
    std::vector<EmPdbResidue> residues;
};

struct EmPdbStrand {
    std::string id;
    int sense;
};
struct EmPdbSheet {
    std::string id;
    std::vector<EmPdbStrand> strands;
};
struct EmPdbHelix {
    std::string id;
    std::string name;
    PdbHelixType type;
    std::string comment;
    int length;
};

struct EmPdbModel {
    std::string id;

    // Primary Structures
    std::vector<EmPdbChain> chains;

    // Secondary Structure
    std::vector<EmPdbHelix> helices;
    std::vector<EmPdbSheet> sheets;
};

struct EmPdbInfo {
    std::vector<EmPdbModel> models;
};
EmPdbInfo emGetPdbInfo() {
    EmPdbInfo info;
    for (auto &p_models : pdb.models) {
        auto &model      = p_models.second;
        auto &model_info = info.models.emplace_back();

        model_info.id = std::to_string(model.serial);
        // Primary Structures
        for (auto &p_chains : model.chains) {
            auto &chain      = p_chains.second;
            auto &chain_info = model_info.chains.emplace_back();

            chain_info.id = std::string(1, chain.chain_id);
            for (auto &residue_ptr : chain.residues) {
                auto &residue_info = chain_info.residues.emplace_back();

                residue_info.id   = std::to_string(residue_ptr->res_seq);
                residue_info.name = std::string(residue_ptr->res_name);
            }
        }

        // Secondary Structure
        for (auto &p_helix : model.helices) {
            auto &helix      = p_helix.second;
            auto &helix_info = model_info.helices.emplace_back();

            helix_info.id      = std::to_string(helix.serial);
            helix_info.name    = std::string(helix.id_code);
            helix_info.comment = std::string(helix.comment);
            helix_info.type    = helix.type;
            helix_info.length  = helix.seq_length;
        }
        for (auto &p_sheet : model.sheets) {
            auto &sheet      = p_sheet.second;
            auto &sheet_info = model_info.sheets.emplace_back();

            sheet_info.id = std::string(sheet.sheet_id);
            for (auto &p_strand : sheet.strands) {
                auto &strand      = p_strand.second;
                auto &strand_info = sheet_info.strands.emplace_back();

                strand_info.id    = std::to_string(strand.strand_id);
                strand_info.sense = strand.sense;
            }
        }
    }

    return info;
}

MolFile emGetMolInfo() {
    return mol;
}

// Binding c++ functions to js api
EMSCRIPTEN_BINDINGS(module) {
    // Define datatype bindings
    emscripten::register_vector<EmPdbResidue>("vector<Residue>");
    emscripten::register_vector<EmPdbChain>("vector<Chain>");
    emscripten::register_vector<EmPdbStrand>("vector<Strand>");
    emscripten::register_vector<EmPdbSheet>("vector<Sheet>");
    emscripten::register_vector<EmPdbHelix>("vector<Helix>");
    emscripten::register_vector<EmPdbModel>("vector<PdbModel>");
    emscripten::register_vector<MolAtom>("vector<Atom>");
    emscripten::register_vector<MolBond>("vector<Bond>");

    // Enums
    emscripten::enum_<Camera::Type>("CameraType")
        .value("PERSPECTIVE", Camera::Type::PERSPECTIVE)
        .value("ORTHOGRAPHIC", Camera::Type::ORTHOGRAPHIC);
    emscripten::enum_<DrawSettings::Mode>("DrawSettingsMode")
        .value("NORMAL", DrawSettings::Mode::NORMAL)
        .value("GOOCH", DrawSettings::Mode::GOOCH)
        .value("BLINN_PHONG", DrawSettings::Mode::BLINN_PHONG);
    emscripten::enum_<MolBondType>("BondType")
        .value("UNKNOWN", MolBondType::UNKNOWN)
        .value("SINGLE", MolBondType::SINGLE)
        .value("DOUBLE", MolBondType::DOUBLE)
        .value("TRIPLE", MolBondType::TRIPLE)
        .value("AROMATIC", MolBondType::AROMATIC)
        .value("SINGLE_OR_DOUBLE", MolBondType::SINGLE_OR_DOUBLE)
        .value("SINGLE_OR_AROMATIC", MolBondType::SINGLE_OR_AROMATIC)
        .value("DOUBLE_OR_AROMATIC", MolBondType::DOUBLE_OR_AROMATIC)
        .value("ANY", MolBondType::ANY);
    emscripten::enum_<PdbResidueColorMode>("ResidueColorMode")
        .value("CHAIN", PdbResidueColorMode::CHAIN)
        .value("SECONDARY", PdbResidueColorMode::SECONDARY)
        .value("AMINO_ACID", PdbResidueColorMode::AMINO_ACID);
    emscripten::enum_<PdbHelixType>("HelixType")
        .value("RIGHT_HANDED_ALPHA", PdbHelixType::RIGHT_HANDED_ALPHA)
        .value("RIGHT_HANDED_OMEGA", PdbHelixType::RIGHT_HANDED_OMEGA)
        .value("RIGHT_HANDED_PI", PdbHelixType::RIGHT_HANDED_PI)
        .value("RIGHT_HANDED_GAMMA", PdbHelixType::RIGHT_HANDED_GAMMA)
        .value("RIGHT_HANDED_310", PdbHelixType::RIGHT_HANDED_310)
        .value("RIGHT_HANDED_ALPHA", PdbHelixType::LEFT_HANDED_ALPHA)
        .value("RIGHT_HANDED_OMEGA", PdbHelixType::LEFT_HANDED_OMEGA)
        .value("RIGHT_HANDED_PI", PdbHelixType::LEFT_HANDED_GAMMA)
        .value("RIGHT_HANDED_GAMMA", PdbHelixType::HELIX_27)
        .value("RIGHT_HANDED_310", PdbHelixType::POLYPROLINE);

    // Structs
    emscripten::value_object<glm::fvec3>("glm::fvec3")
        .field("x", &glm::fvec3::x)
        .field("y", &glm::fvec3::y)
        .field("z", &glm::fvec3::z);
    emscripten::value_object<glm::fvec4>("glm::fvec4")
        .field("x", &glm::fvec4::x)
        .field("y", &glm::fvec4::y)
        .field("z", &glm::fvec4::z)
        .field("w", &glm::fvec4::w);
    // @todo seperate bond and atom type
    emscripten::value_object<MolAtom>("Atom")
        .field("position", &MolAtom::position)
        .field("symbol", &MolAtom::symbol)
        .field("massDifference", &MolAtom::mass_difference)
        .field("charge", &MolAtom::charge)
        .field("hydrogenCount", &MolAtom::hydrogen_count)
        .field("valence", &MolAtom::valence);
    emscripten::value_object<MolBond>("Bond")
        .field("atom1Id", &MolBond::atom_1_index)
        .field("atom2Id", &MolBond::atom_2_index)
        .field("type", &MolBond::type);
    emscripten::value_object<MolFile>("MolInfo")
        .field("title", &MolFile::title)
        .field("comments", &MolFile::comments)
        .field("metadata", &MolFile::metadata)
        .field("atoms", &MolFile::atoms)
        .field("bonds", &MolFile::bonds);
    emscripten::value_object<EmPdbResidue>("Residue")
        .field("name", &EmPdbResidue::name)
        .field("id", &EmPdbResidue::id);
    emscripten::value_object<EmPdbChain>("Chain")
        .field("residues", &EmPdbChain::residues)
        .field("id", &EmPdbChain::id);
    emscripten::value_object<EmPdbStrand>("Strand")
        .field("sense", &EmPdbStrand::sense)
        .field("id", &EmPdbStrand::id);
    emscripten::value_object<EmPdbSheet>("Sheet")
        .field("strands", &EmPdbSheet::strands)
        .field("id", &EmPdbSheet::id);
    emscripten::value_object<EmPdbHelix>("Helix")
        .field("type", &EmPdbHelix::type)
        .field("name", &EmPdbHelix::name)
        .field("comment", &EmPdbHelix::comment)
        .field("length", &EmPdbHelix::length)
        .field("id", &EmPdbHelix::id);
    emscripten::value_object<EmPdbModel>("PdbModel")
        .field("sheets", &EmPdbModel::sheets)
        .field("helices", &EmPdbModel::helices)
        .field("chains", &EmPdbModel::chains)
        .field("id", &EmPdbModel::id);
    emscripten::value_object<EmPdbInfo>("PdbInfo")
        .field("models", &EmPdbInfo::models);
    emscripten::value_object<DrawSettings>("DrawSettings")
        .field("clearColor", &DrawSettings::clear_color)
        .field("lightColor", &DrawSettings::light_color)
        .field("lightDirection", &DrawSettings::light_direction)
        .field("mode", &DrawSettings::mode);
    emscripten::value_object<PdbDrawSettings>("PdbDrawSettings")
        .field("drawHeteroAtoms", &PdbDrawSettings::draw_hetero_atoms)
        .field("drawWaterAtoms", &PdbDrawSettings::draw_water_atoms)
        .field("drawResidueAtoms", &PdbDrawSettings::draw_residue_atoms)
        .field("drawResidueRibbons", &PdbDrawSettings::draw_residue_ribbons)
        .field("heteroAtomsAlpha", &PdbDrawSettings::hetero_atoms_alpha)
        .field("waterAtomsAlpha", &PdbDrawSettings::water_atoms_alpha)
        .field("residueAtomsAlpha", &PdbDrawSettings::residue_atoms_alpha)
        .field("residueRibbonsAlpha", &PdbDrawSettings::residue_ribbons_alpha)
        .field("residueColorMode", &PdbDrawSettings::residue_color_mode);

    // Define functional bindings
    emscripten::function("loadPdb", &emLoadPdb);
    emscripten::function("loadMol", &emLoadMol);
    emscripten::function("loadPdbDictionary", &emLoadPdbDictionary);

    emscripten::function("getDrawSettings", &emGetDrawSettings);
    emscripten::function("getPdbDrawSettings", &emGetPdbDrawSettings);
    emscripten::function("getCameraType", &emGetCameraType);
    emscripten::function("getPdbInfo", &emGetPdbInfo);
    emscripten::function("getMolInfo", &emGetMolInfo);

    emscripten::function("setDrawSettings", &emSetDrawSettings);
    emscripten::function("setPdbDrawSettings", &emSetPdbDrawSettings);
    emscripten::function("setCameraType", &emSetCameraType);
    emscripten::function("setViewportSize", &emSetViewportSize);
}

EM_BOOL emScrollCallback(int eventType, const EmscriptenWheelEvent *wheelEvent, void *userData) {
    // std::cout << "Mousemove event, x: " << mouseEvent->clientX << ", y: " <<  mouseEvent->clientY  << "\n";
    double xoffset = wheelEvent->deltaX / 100.0f;
    double yoffset = wheelEvent->deltaY / 100.0f;
    if (Controls::scroll_offset.x != xoffset || Controls::scroll_offset.y != yoffset)
        Controls::scrolled = true;
    else
        Controls::scrolled = false;

    Controls::scroll_offset.x = xoffset;
    Controls::scroll_offset.y = yoffset;
    return EM_TRUE;
}

EM_BOOL emMouseDownCallback(int eventType, const EmscriptenMouseEvent *mouseEvent, void *userData) {
    if (mouseEvent->button == 0)
        Controls::left_mouse = true;
    else if (mouseEvent->button == 2)
        Controls::right_mouse = true;
    return EM_TRUE;
}
EM_BOOL emMouseMoveCallback(int eventType, const EmscriptenMouseEvent *mouseEvent, void *userData) {
    Controls::old_mouse_position = Controls::mouse_position;
    Controls::mouse_position     = glm::vec2(mouseEvent->targetX, mouseEvent->targetY);
    Controls::mouse_moving       = true;
    return EM_TRUE;
}
EM_BOOL emMouseUpCallback(int eventType, const EmscriptenMouseEvent *mouseEvent, void *userData) {
    if (mouseEvent->button == 0)
        Controls::left_mouse = false;
    else if (mouseEvent->button == 2)
        Controls::right_mouse = false;
    return EM_TRUE;
}

EM_BOOL emTouchStartCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData) {
    Controls::touched = touchEvent->numTouches == 1;
    if (Controls::touched) {
        Controls::mouse_position     = glm::vec2(touchEvent->touches[0].targetX, touchEvent->touches[0].targetY);
        Controls::old_mouse_position = Controls::mouse_position;
    }
    // Allow focus to return
    return EM_FALSE;
}
EM_BOOL emTouchMoveCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData) {
    if (Controls::touched) {
        Controls::old_mouse_position = Controls::mouse_position;
        Controls::mouse_position     = glm::dvec2(touchEvent->touches[0].targetX, touchEvent->touches[0].targetY);
    }

    Controls::touched = touchEvent->numTouches == 1;
    return EM_FALSE;
}
EM_BOOL emTouchEndCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData) {
    Controls::touched = false;
    // Allow focus to return
    return EM_FALSE;
}
EM_BOOL emTouchCancelCallback(int eventType, const EmscriptenTouchEvent *touchEvent, void *userData) {
    Controls::touched = false;
    // Allow focus to return
    return EM_FALSE;
}

#else

void glfwScrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    if (Controls::scroll_offset.x != xoffset || Controls::scroll_offset.y != yoffset)
        Controls::scrolled = true;
    else
        Controls::scrolled = false;

    Controls::scroll_offset.x = xoffset;
    Controls::scroll_offset.y = yoffset;
}

void glfwCursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    glm::dvec2 new_mouse_position(xpos, ypos);
    if (Controls::mouse_position.x < 0)
        Controls::mouse_position = new_mouse_position;
    Controls::old_mouse_position = Controls::mouse_position;
    Controls::mouse_position     = new_mouse_position;

    Controls::mouse_moving = true;
}

void glfwResizeCallback(GLFWwindow* window, int width, int height) {
    window_width  = width;
    window_height = height;
    camera.set_aspect_ratio((double)window_width / (double)window_height);
    glfwSetWindowSize(window, window_width, window_height);
}

void glfwMouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    Controls::right_mouse = button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS;
    Controls::left_mouse  = button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS;
}

#endif

[[nodiscard]] bool init() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLEW\n";
        return false;
    }

    // Create main window.
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT,
                   GL_TRUE);    // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, 1);

#if EMSCRIPTEN
    glfwWindowHint(GLFW_MAXIMIZED, 1);
#else
    glfwWindowHint(GLFW_RESIZABLE, 1);

    window_width  = 1024;
    window_height = 700;
#endif

    window = glfwCreateWindow(window_width, window_height, "cviz", NULL, NULL);
    if (window == NULL) {
        std::cerr << "Failed to create glfw window\n";
        return false;
    }
    glfwMakeContextCurrent(window);

#ifndef EMSCRIPTEN
    //
    // Init OpenGL
    //
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW\n";
        return false;
    }
    glEnable(GL_MULTISAMPLE);
#endif

    //
    // Setup callbacks
    //
#if EMSCRIPTEN
    emscripten_set_wheel_callback("#canvas", 0, true, emScrollCallback);
    emscripten_set_mousedown_callback("#canvas", 0, true, emMouseDownCallback);
    emscripten_set_mousemove_callback("#canvas", 0, true, emMouseMoveCallback);
    emscripten_set_mouseup_callback("#canvas", 0, true, emMouseUpCallback);

    emscripten_set_touchstart_callback("#canvas", 0, true, emTouchStartCallback);
    emscripten_set_touchmove_callback("#canvas", 0, true, emTouchMoveCallback);
    emscripten_set_touchend_callback("#canvas", 0, true, emTouchEndCallback);
    emscripten_set_touchcancel_callback("#canvas", 0, true, emTouchCancelCallback);
#else
    glfwSetWindowSizeCallback(window, &glfwResizeCallback);
    glfwSetScrollCallback(window, &glfwScrollCallback);
    glfwSetMouseButtonCallback(window, &glfwMouseButtonCallback);
    glfwSetCursorPosCallback(window, &glfwCursorPosCallback);

    camera.set_aspect_ratio((double)window_width / (double)window_height);
#endif

    return init_graphics();
}
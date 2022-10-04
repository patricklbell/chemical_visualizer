#ifndef EDITOR_H
#define EDITOR_H

#include <vector>

#include "entities.hpp"
#include "controls.hpp"
#include "loader.hpp"

void loadPdbFileHelper(std::string path, Entities &entities, Camera &camera);
void loadMolFileHelper(std::string path, Entities &entities, Camera &camera);

void initGui();
void drawGui(Camera &camera, Entities &entities);

enum class UiMode {
	NONE = 0,
	PDB,
	MOL,
};

namespace ui {
#if defined(EMSCRIPTEN)
    extern bool fullscreen;
#endif
    extern UiMode mode;
	extern MolFile molfile;

	extern PdbDictionary pdb_dictionary;
	extern PdbFile pdbfile;
	extern PdbDrawSettings pdbfile_settings;
    extern int selected_pdb_model_serial;

	extern std::string loaded_file_path;
	extern bool dark_mode;
}

#endif

#ifndef EDITOR_H
#define EDITOR_H

#include <vector>

#include "entities.hpp"
#include "controls.hpp"
#include "loader.hpp"

void initGui();
void drawGui(Camera &camera, Entities &entities);

enum class UiMode {
	NONE = 0,
	PDB,
	MOL,
};

namespace ui {
    extern UiMode mode;
	extern MolFile molfile;

	extern PdbFile pdbfile;
	extern PdbDrawSettings pdbfile_settings;

	extern std::string loaded_file_path;
}

#endif

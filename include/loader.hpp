#ifndef LOADER_HPP
#define LOADER_HPP

#include <vector>

#include "entities.hpp"

// http://c4.cabrillo.edu/404/ctfile.pdf
struct MolAtom {
    glm::fvec3 position;
    char symbol[3];
    // Difference of this atom from mass in periodic table
    float mass_difference = 0;
    float charge = 0;
    // n + 1 -> Hn meaning n extra hydrogen atoms are allowed to be drawn
    int hydrogen_count = 1;
    // Number of bonds to this atom
    int valence = 0;
};

enum class MolBondType : unsigned int {
    UNKNOWN             = 0,
    SINGLE              = 1,
    DOUBLE              = 2,
    TRIPLE              = 3,
    AROMATIC            = 4,
    SINGLE_OR_DOUBLE    = 5,
    SINGLE_OR_AROMATIC  = 6,
    DOUBLE_OR_AROMATIC  = 7,
    ANY                 = 8,
};

//enum class MolBondStereo : unsigned int {
//    NOT     = 1,
//    UP      = 2,
//    DOWN    = 3,
//    EITHER  = 4,
//    EITHER_DOUBLE = 5,
//};
//enum class MolBondTopology : unsigned int {
//    EITHER  = 0,
//    RING    = 1,
//    CHAIN   = 2,
//};

struct MolBond {
    int atom_1_index;
    int atom_2_index;
    MolBondType type;
};

struct MolFile {
    char title[80];
    char comments[80];
    int num_atoms;
    MolAtom *atoms;
    int num_bonds;
    MolBond *bonds;
};

void loadMolfile(MolFile &data, std::string path);
glm::vec3 createEntitiesFromMolfile(std::vector<Entity*> &entities, MolFile &data);

#endif

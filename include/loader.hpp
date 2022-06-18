#ifndef LOADER_HPP
#define LOADER_HPP

#include <vector>
#include <unordered_map>

#include "entities.hpp"

// http://c4.cabrillo.edu/404/ctfile.pdf
struct MolAtom {
    glm::fvec3 position;
    char symbol[3] = "";
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

void loadMolFile(MolFile &data, std::string path);
glm::vec3 createEntitiesFromMolFile(std::vector<Entity*> &entities, MolFile &data);

struct PdbAtom {
    int serial;
    char symbol[3] = "";
    //char name[5];
    //char alternate_location;
    //char residue_name[5];
    //int residue_sequence;
    //char chain_id;
    glm::fvec3 position; // Orthogonal coordinates in Angstroms
};

enum class PdbBondType : unsigned int {
    UNKNOWN             = 0,
    SINGLE              = 1,
    DOUBLE              = 2,
    TRIPLE              = 3,
};

struct PdbBond {
    int atom_1_index;
    int atom_2_index;
    PdbBondType type;
};

struct PdbModel {
    int serial = 0;
    std::vector<PdbAtom> atoms;
    std::vector<PdbBond> bonds;
};

struct PdbFile {
    PdbModel heterogen_model;

    std::vector<PdbModel> polymer_models;
};

void loadPdbFile(PdbFile &data, std::string path);
glm::vec3 createEntitiesFromPdbFile(std::vector<Entity*> &entities, PdbFile &data);

// --------------------------------Color LUTs-------------------------------- //

const std::unordered_map<std::string, glm::vec3> symbol_to_color_lut = {
    {"",   glm::vec3(255,255,255)},
    {"H",  glm::vec3(55,255,255)},
    {"He", glm::vec3(217,255,255)},
    {"Li", glm::vec3(204,128,255)},
    {"Be", glm::vec3(194,255,0)},
    {"B",  glm::vec3(55,181,181)},
    {"C",  glm::vec3(44,144,144)},
    {"N",  glm::vec3(8,80,248)},
    {"O",  glm::vec3(55,13,13)},
    {"F",  glm::vec3(44,224,80)},
    {"Ne", glm::vec3(179,227,245)},
    {"Na", glm::vec3(171,92,242)},
    {"Mg", glm::vec3(138,255,0)},
    {"Al", glm::vec3(191,166,166)},
    {"Si", glm::vec3(240,200,160)},
    {"P",  glm::vec3(55,128,0)},
    {"S",  glm::vec3(55,255,48)},
    {"Cl", glm::vec3(31,240,31)},
    {"Ar", glm::vec3(128,209,227)},
    {"K",  glm::vec3(43,64,212)},
    {"Ca", glm::vec3(61,255,0)},
    {"Sc", glm::vec3(230,230,230)},
    {"Ti", glm::vec3(191,194,199)},
    {"V",  glm::vec3(66,166,171)},
    {"Cr", glm::vec3(138,153,199)},
    {"Mn", glm::vec3(156,122,199)},
    {"Fe", glm::vec3(224,102,51)},
    {"Co", glm::vec3(240,144,160)},
    {"Ni", glm::vec3(80,208,80)},
    {"Cu", glm::vec3(200,128,51)},
    {"Zn", glm::vec3(125,128,176)},
    {"Ga", glm::vec3(194,143,143)},
    {"Ge", glm::vec3(102,143,143)},
    {"As", glm::vec3(189,128,227)},
    {"Se", glm::vec3(255,161,0)},
    {"Br", glm::vec3(166,41,41)},
    {"Kr", glm::vec3(92,184,209)},
    {"Rb", glm::vec3(112,46,176)},
    {"Sr", glm::vec3(0,255,0)},
    {"Y",  glm::vec3(48,255,255)},
    {"Zr", glm::vec3(148,224,224)},
    {"Nb", glm::vec3(115,194,201)},
    {"Mo", glm::vec3(84,181,181)},
    {"Tc", glm::vec3(59,158,158)},
    {"Ru", glm::vec3(36,143,143)},
    {"Rh", glm::vec3(10,125,140)},
    {"Pd", glm::vec3(0,105,133)},
    {"Ag", glm::vec3(192,192,192)},
    {"Cd", glm::vec3(255,217,143)},
    {"In", glm::vec3(166,117,115)},
    {"Sn", glm::vec3(102,128,128)},
    {"Sb", glm::vec3(158,99,181)},
    {"Te", glm::vec3(212,122,0)},
    {"I",  glm::vec3(48,0,148)},
    {"Xe", glm::vec3(66,158,176)},
    {"Cs", glm::vec3(87,23,143)},
    {"Ba", glm::vec3(0,201,0)},
    {"La", glm::vec3(112,212,255)},
    {"Ce", glm::vec3(255,255,199)},
    {"Pr", glm::vec3(217,255,199)},
    {"Nd", glm::vec3(199,255,199)},
    {"Pm", glm::vec3(163,255,199)},
    {"Sm", glm::vec3(143,255,199)},
    {"Eu", glm::vec3(97,255,199)},
    {"Gd", glm::vec3(69,255,199)},
    {"Tb", glm::vec3(48,255,199)},
    {"Dy", glm::vec3(31,255,199)},
    {"Ho", glm::vec3(0,255,156)},
    {"Er", glm::vec3(0,230,117)},
    {"Tm", glm::vec3(0,212,82)},
    {"Yb", glm::vec3(0,191,56)},
    {"Lu", glm::vec3(0,171,36)},
    {"Hf", glm::vec3(77,194,255)},
    {"Ta", glm::vec3(77,166,255)},
    {"W",  glm::vec3(3,148,214)},
    {"Re", glm::vec3(38,125,171)},
    {"Os", glm::vec3(38,102,150)},
    {"Ir", glm::vec3(23,84,135)},
    {"Pt", glm::vec3(208,208,224)},
    {"Au", glm::vec3(255,209,35)},
    {"Hg", glm::vec3(184,184,208)},
    {"Tl", glm::vec3(166,84,77)},
    {"Pb", glm::vec3(87,89,97)},
    {"Bi", glm::vec3(158,79,181)},
    {"Po", glm::vec3(171,92,0)},
    {"At", glm::vec3(117,79,69)},
    {"Rn", glm::vec3(66,130,150)},
    {"Fr", glm::vec3(66,0,102)},
    {"Ra", glm::vec3(0,125,0)},
    {"Ac", glm::vec3(112,171,250)},
    {"Th", glm::vec3(0,186,255)},
    {"Pa", glm::vec3(0,161,255)},
    {"U",  glm::vec3(0,143,255)},
    {"Np", glm::vec3(0,128,255)},
    {"Pu", glm::vec3(0,107,255)},
    {"Am", glm::vec3(84,92,242)},
    {"Cm", glm::vec3(120,92,227)},
    {"Bk", glm::vec3(138,79,227)},
    {"Cf", glm::vec3(161,54,212)},
    {"Es", glm::vec3(179,31,212)},
    {"Fm", glm::vec3(179,31,186)},
    {"Md", glm::vec3(179,13,166)}
};

#endif

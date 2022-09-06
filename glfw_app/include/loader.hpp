#ifndef LOADER_HPP
#define LOADER_HPP

#include <vector>
#include <unordered_map>

#include "entities.hpp"
#include "graphics.hpp"

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
void createEntitiesFromMolFile(Entities &entities, MolFile &data, Camera &camera);

struct PdbAtom {
    // -1 indicates heterogen_model, otherwise index of model
    int model_num;
    bool is_heterogen;

    int serial;
    char name[5] = "";

    int  res_seq;
    char res_name[4];
    char chain_id;
    char i_code;

    glm::fvec3 position; // Orthogonal coordinates in Angstroms
    
    char symbol[3] = "";
};

enum class PdbConnectionType : unsigned int {
    UNKNOWN             = 0,
    SINGLE              = 1,
    SINGLE_REDUNDANT    = 2,
    DOUBLE              = 3,
    DOUBLE_REDUNDANT    = 4,
    TRIPLE              = 5,
    TRIPLE_REDUNDANT    = 6,
};

struct PdbConnection {
    int atom_1_id;
    int atom_2_id;
    PdbConnectionType type;
};

enum class PdbHelixType : unsigned int {
    RIGHT_HANDED_ALPHA  = 1,
    RIGHT_HANDED_OMEGA  = 2,
    RIGHT_HANDED_PI     = 3,
    RIGHT_HANDED_GAMMA  = 4,
    RIGHT_HANDED_310    = 5,
    LEFT_HANDED_ALPHA   = 6,
    LEFT_HANDED_OMEGA   = 7,
    LEFT_HANDED_GAMMA   = 8,
    HELIX_27            = 9,
    POLYPROLINE         = 10,
};

enum class PdbResidueType : unsigned int {
    UNKNOWN = 0,
    COIL    = 1,
    HELIX   = 2,
    STRAND  = 3,
};

struct PdbHelix {
    int serial;
    char id_code[4];

    int  init_seq_num;
    char init_res_name[4];
    char init_chain_id;
    char init_i_code;

    int  end_seq_num;
    char end_res_name[4];
    char end_chain_id;
    char end_i_code;

    int seq_length;

    PdbHelixType type;
};

struct PdbStrand {
    int strand_id;
    int sense;

    int  init_seq_num;
    char init_res_name[4];
    char init_chain_id;
    char init_i_code;

    int  end_seq_num;
    char end_res_name[4];
    char end_chain_id;
    char end_i_code;

    // If first strand cur and prev ignored
    bool is_first;

    char cur_atom_name[5];
    int  cur_seq_num;
    char cur_res_name[4];
    char cur_chain_id;
    char cur_i_code;

    char prev_atom_name[5];
    int  prev_seq_num;
    char prev_res_name[4];
    char prev_chain_id;
    char prev_i_code;
};

struct PdbSheet {
    char sheet_id[4];
    std::unordered_map<int, PdbStrand> strands; // strand_id -> strand

    int num_strands;
};

struct PdbResidue {
    int  res_seq;
    char res_name[4];
    char chain_id;
    char i_code;

    std::vector<int> atom_ids;
    // @note LINK should guaranteed that atom_name is unique in residue
    std::unordered_map<std::string, int> atom_name_id; // atom_name -> serial

    // Used for constructing secondary structures, not inherent to residue
    // assume coil as pdb doesn't specify
    PdbResidueType type = PdbResidueType::COIL;
};

// A contiguous sequence of residues (in secondary structures) which form one mesh
struct PdbChain {
    char chain_id;
    std::vector<PdbResidue *> residues;
};

struct PdbModel {
    int serial = 0;

    // Primary structures
    std::unordered_map<int,  PdbAtom> atoms;  // serial -> atom
    std::unordered_map<long, PdbConnection> connections; // hash(serial1, serial2) -> bond

    // @note c++ doesn't seem to pre define char hash function
    std::unordered_map<int, PdbChain> chains; // chain_id -> chain
    
    // Secondary structures
    std::unordered_map<int, PdbHelix> helices; // serial -> helix // @note this may cause collisions?
    std::unordered_map<std::string, PdbSheet> sheets; // sheet_id -> sheet
    
    std::unordered_map<int, PdbResidue> residues; // res_seq -> residue
};

struct PdbFile {
    // Each model represents the same structure, mainly for NMR entries
    std::vector<PdbModel> models;
};

// Approximates plane of residue with peptide bonds
struct PeptidePlane {
    PdbResidue *residue_1;
    PdbResidue *residue_2;

    glm::vec3 position, normal, forward, right;

    // Peptide bond could be flipped relative to previous residue
    bool flipped = false;
};

void loadPdbFile(PdbFile &data, std::string path);
void createEntitiesFromPdbFile(Entities &entities, PdbFile &data, Camera &camera);

// --------------------------------Color LUTs-------------------------------- //

const std::unordered_map<std::string, glm::vec3> symbol_to_color_lut = {
    {"",   glm::vec3(255,255,255)},
    {"H",  glm::vec3(55,255,255)},
    {"He", glm::vec3(217,255,255)},
    {"Li", glm::vec3(204,128,255)},
    {"Be", glm::vec3(194,255,0)},
    {"B",  glm::vec3(55,181,181)},
    {"C",  glm::vec3(100,100,100)},
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

#ifndef LOADER_HPP
#define LOADER_HPP

#include <vector>
#include <unordered_map>

#include "entities.hpp"
#include "controls.hpp"

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
    
    std::vector<MolAtom> atoms;
    std::vector<MolBond> bonds;

    void clear();
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

    // 
    // Used for rendering 
    //
    int entity_id = -1;
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

    // 
    // Used for rendering 
    //
    int entity_id = -1; // Id of last entity
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
    // assume coil if pdb doesn"t specify
    PdbResidueType type = PdbResidueType::COIL;

    // 
    // Used for rendering 
    //
    int entity_id = -1;
    bool mesh_generated = false;
    Mesh mesh; // The first residue of the chain owns the mesh
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

    // @note c++ doesn"t seem to pre define char hash function
    std::unordered_map<int, PdbChain> chains; // chain_id -> chain
    
    // Secondary structures
    std::unordered_map<int, PdbHelix> helices; // serial -> helix // @note this may cause collisions?
    std::unordered_map<std::string, PdbSheet> sheets; // sheet_id -> sheet
    
    std::unordered_map<int, PdbResidue> residues; // res_seq -> residue

    bool renderable = false; // Whether meshes have been generated
};

struct PdbFile {
    // Each model represents the same structure, mainly for NMR entries
    std::unordered_map<int, PdbModel> models; // serial -> model

    void clear();
};

struct PdbDictionaryConnect {
    std::string atom_name_1 = "";
    std::string atom_name_2 = "";
    PdbConnectionType type = PdbConnectionType::UNKNOWN;
};

struct PdbDictionary {
    // This is a map between residue names and the residue"s connections, (stored as map between atom_name_1+2) -> connection
    // We need to use names since the dictionary doesn"t have a concept of atom id
    std::unordered_map<std::string, std::unordered_map<std::string, PdbDictionaryConnect>> residues; 
};

enum PdbResidueColorMode : int {
    CHAIN = 0,
    SECONDARY,
    AMINO_ACID,
    NUM_MODES,
};

struct PdbDrawSettings {
    bool draw_hetero_atoms    = true;
    bool draw_water_atoms     = false;
    bool draw_residue_atoms   = false;
    bool draw_residue_ribbons = true;

    // Alpha channel of color
    float hetero_atoms_alpha    = 1.0;
    float water_atoms_alpha     = 1.0;
    float residue_atoms_alpha   = 1.0;
    float residue_ribbons_alpha = 1.0;

    PdbResidueColorMode residue_color_mode;
};

// Approximates plane of residue with peptide bonds
struct PeptidePlane {
    PdbResidue *residue_1;
    PdbResidue *residue_2;

    glm::vec3 position, normal, forward, right;
};

void createDebugCartesian(const glm::vec3 &p, const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c, Entities &entities, float r);

void loadPdbDictionaryFile(PdbDictionary &dict, std::string_view path);
void loadPdbFile(PdbFile &data, std::string path, PdbDictionary *dict=nullptr);
void updateEntityColorsFromPdbModel(Entities& entities, PdbModel& model, PdbDrawSettings& settings);
void createPdbModelMeshes(PdbModel& model);
void createEntitiesFromPdbModel(Entities& entities, PdbModel& model, PdbDrawSettings& settings, Camera& camera);

// --------------------------------Color LUTs-------------------------------- //
// Based on: https://jmol.sourceforge.net/jscolors/

const std::unordered_map<std::string, glm::vec3> symbol_to_color_lut = {
    {"",   glm::vec3(255,255,255)},
    {"H",  glm::vec3(55,255,255)},
    {"HE", glm::vec3(217,255,255)},
    {"LI", glm::vec3(204,128,255)},
    {"BE", glm::vec3(194,255,0)},
    {"B",  glm::vec3(55,181,181)},
    {"C",  glm::vec3(100,100,100)},
    {"N",  glm::vec3(8,80,248)},
    {"O",  glm::vec3(55,13,13)},
    {"F",  glm::vec3(44,224,80)},
    {"NE", glm::vec3(179,227,245)},
    {"NA", glm::vec3(171,92,242)},
    {"MG", glm::vec3(138,255,0)},
    {"AL", glm::vec3(191,166,166)},
    {"SI", glm::vec3(240,200,160)},
    {"P",  glm::vec3(55,128,0)},
    {"S",  glm::vec3(55,255,48)},
    {"CL", glm::vec3(31,240,31)},
    {"AR", glm::vec3(128,209,227)},
    {"K",  glm::vec3(43,64,212)},
    {"CA", glm::vec3(61,255,0)},
    {"SC", glm::vec3(230,230,230)},
    {"TI", glm::vec3(191,194,199)},
    {"V",  glm::vec3(66,166,171)},
    {"CR", glm::vec3(138,153,199)},
    {"MN", glm::vec3(156,122,199)},
    {"FE", glm::vec3(224,102,51)},
    {"CO", glm::vec3(240,144,160)},
    {"NI", glm::vec3(80,208,80)},
    {"CU", glm::vec3(200,128,51)},
    {"ZN", glm::vec3(125,128,176)},
    {"GA", glm::vec3(194,143,143)},
    {"GE", glm::vec3(102,143,143)},
    {"AS", glm::vec3(189,128,227)},
    {"SE", glm::vec3(255,161,0)},
    {"BR", glm::vec3(166,41,41)},
    {"KR", glm::vec3(92,184,209)},
    {"RB", glm::vec3(112,46,176)},
    {"SR", glm::vec3(0,255,0)},
    {"Y",  glm::vec3(48,255,255)},
    {"ZR", glm::vec3(148,224,224)},
    {"NB", glm::vec3(115,194,201)},
    {"MO", glm::vec3(84,181,181)},
    {"TC", glm::vec3(59,158,158)},
    {"RU", glm::vec3(36,143,143)},
    {"RH", glm::vec3(10,125,140)},
    {"PD", glm::vec3(0,105,133)},
    {"AG", glm::vec3(192,192,192)},
    {"CD", glm::vec3(255,217,143)},
    {"IN", glm::vec3(166,117,115)},
    {"SN", glm::vec3(102,128,128)},
    {"SB", glm::vec3(158,99,181)},
    {"TE", glm::vec3(212,122,0)},
    {"I",  glm::vec3(48,0,148)},
    {"XE", glm::vec3(66,158,176)},
    {"CS", glm::vec3(87,23,143)},
    {"BA", glm::vec3(0,201,0)},
    {"LA", glm::vec3(112,212,255)},
    {"CE", glm::vec3(255,255,199)},
    {"PR", glm::vec3(217,255,199)},
    {"ND", glm::vec3(199,255,199)},
    {"PM", glm::vec3(163,255,199)},
    {"SM", glm::vec3(143,255,199)},
    {"EU", glm::vec3(97,255,199)},
    {"GD", glm::vec3(69,255,199)},
    {"TB", glm::vec3(48,255,199)},
    {"DY", glm::vec3(31,255,199)},
    {"HO", glm::vec3(0,255,156)},
    {"ER", glm::vec3(0,230,117)},
    {"TM", glm::vec3(0,212,82)},
    {"YB", glm::vec3(0,191,56)},
    {"LU", glm::vec3(0,171,36)},
    {"HF", glm::vec3(77,194,255)},
    {"TA", glm::vec3(77,166,255)},
    {"W",  glm::vec3(3,148,214)},
    {"RE", glm::vec3(38,125,171)},
    {"OS", glm::vec3(38,102,150)},
    {"IR", glm::vec3(23,84,135)},
    {"PT", glm::vec3(208,208,224)},
    {"AU", glm::vec3(255,209,35)},
    {"HG", glm::vec3(184,184,208)},
    {"TL", glm::vec3(166,84,77)},
    {"PB", glm::vec3(87,89,97)},
    {"BI", glm::vec3(158,79,181)},
    {"PO", glm::vec3(171,92,0)},
    {"AT", glm::vec3(117,79,69)},
    {"RN", glm::vec3(66,130,150)},
    {"FR", glm::vec3(66,0,102)},
    {"RA", glm::vec3(0,125,0)},
    {"AC", glm::vec3(112,171,250)},
    {"TH", glm::vec3(0,186,255)},
    {"PA", glm::vec3(0,161,255)},
    {"U",  glm::vec3(0,143,255)},
    {"NP", glm::vec3(0,128,255)},
    {"PU", glm::vec3(0,107,255)},
    {"AM", glm::vec3(84,92,242)},
    {"CM", glm::vec3(120,92,227)},
    {"BK", glm::vec3(138,79,227)},
    {"CF", glm::vec3(161,54,212)},
    {"ES", glm::vec3(179,31,212)},
    {"FM", glm::vec3(179,31,186)},
    {"MD", glm::vec3(179,13,166)}
};

const glm::vec3 unknown_residue_color = glm::vec3(190/255.f,160/255.f,110/255.f);
const std::unordered_map<std::string, glm::vec3> residue_to_color_lut = {
    {"ALA",	glm::vec3(200/255.f,200/255.f,200/255.f)},
    {"ARG",	glm::vec3(20/255.f,90/255.f,255  /255.f)},
    {"ASN",	glm::vec3(0/255.f,220/255.f,220  /255.f)},
    {"ASP",	glm::vec3(230/255.f,10/255.f,10  /255.f)},
    {"CYS",	glm::vec3(230/255.f,230/255.f,0  /255.f)},
    {"GLN",	glm::vec3(0/255.f,220/255.f,220  /255.f)},
    {"GLU",	glm::vec3(230/255.f,10/255.f,10  /255.f)},
    {"GLY",	glm::vec3(235/255.f,235/255.f,235/255.f)},
    {"HIS",	glm::vec3(130/255.f,130/255.f,210/255.f)},
    {"ILE",	glm::vec3(15/255.f,130/255.f,15  /255.f)},
    {"LEU",	glm::vec3(15/255.f,130/255.f,15  /255.f)},
    {"LYS",	glm::vec3(20/255.f,90/255.f,255  /255.f)},
    {"MET",	glm::vec3(230/255.f,230/255.f,0  /255.f)},
    {"PHE",	glm::vec3(50/255.f,50/255.f,170  /255.f)},
    {"PRO",	glm::vec3(220/255.f,150/255.f,130/255.f)},
    {"SER",	glm::vec3(250/255.f,150/255.f,0  /255.f)},
    {"THR",	glm::vec3(250/255.f,150/255.f,0  /255.f)},
    {"TRP",	glm::vec3(180/255.f,90/255.f,180 /255.f)},
    {"TYR",	glm::vec3(50/255.f,50/255.f,170  /255.f)},
    {"VAL",	glm::vec3(15/255.f,130/255.f,15  /255.f)},
    {"ASX",	glm::vec3(255/255.f,105/255.f,180/255.f)},
    {"GLX",	glm::vec3(255/255.f,105/255.f,180/255.f)},
};

const auto helix_color = glm::vec4(glm::vec3(255, 0, 128) / 255.f, 1.0);
const auto strand_color = glm::vec4(glm::vec3(255, 200, 0) / 255.f, 1.0);
const auto coil_color = glm::vec4(glm::vec3(96, 128, 255) / 255.f, 1.0);

#endif

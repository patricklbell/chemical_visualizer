#include <unordered_map>

#include <glm/glm.hpp>
#include "assets.hpp"
#include "glm/gtc/quaternion.hpp"

#include "loader.hpp"
#include "entities.hpp"
#include "graphics.hpp"
#include "utilities.hpp"

std::unordered_map<std::string, glm::vec3> symbol_to_color_lut = {
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

void loadMolfile(MolFile &data, std::string path){
    FILE *f;
    f=fopen(path.c_str(), "r");
    if (!f) {
        fprintf(stderr, "Error in reading mesh file %s.\n", path.c_str());
        return;
    }

    printf("----------------Loading Molecule File %s----------------\n", path.c_str());

    fgets(data.title, 80, f);
    printf("Title: %s\n", data.title);
    // @todo read other metadata
    char metadata[80];
    fgets(metadata, 80, f);
    fgets(data.comments, 80, f);
    printf("Comments: %s\n", data.comments);

    fscanf(f, "%d %d %*[^\n] ", &data.num_atoms, &data.num_bonds);
    printf("Num of atoms is %d and bonds is %d\n", data.num_atoms, data.num_bonds);

    data.atoms = (MolAtom*)malloc(sizeof(MolAtom)*data.num_atoms);
    for(int i = 0; i < data.num_atoms; ++i){
        auto &a = data.atoms[i];

        fscanf(f, "%f %f %f", &a.position.x, &a.position.y, &a.position.z);
        fscanf(f, "%s", a.symbol);
        int md;
        fscanf(f, "%d %*d %d %*d %d %*[^\n]", &md, &a.hydrogen_count, &a.valence);
        a.mass_difference = (float)md;

        printf("Atom %s position is x: %f y: %f z: %f\n", a.symbol, a.position.x, a.position.y, a.position.z);
    }

    data.bonds = (MolBond*)malloc(sizeof(MolBond)*data.num_bonds);
    for(int i = 0; i < data.num_bonds; ++i){
        auto &b = data.bonds[i];

        fscanf(f, "%d %d %u %*[^\n]", &b.atom_1_index, &b.atom_2_index, &b.type);
        b.atom_1_index--;
        b.atom_2_index--;
        printf("Bond %d from atom %d --> %d\n", i, b.atom_1_index, b.atom_2_index);
    }
}

glm::vec3 createEntitiesFromMolfile(std::vector<Entity*> &entities, MolFile &data){
    static const float sphere_r = 0.2;
    static const float cylinder_r = 0.035;
    static const float bond_gap_r = 0.08;
    //static const float cylinder_offset = std::sqrt(sphere_r*sphere_r - cylinder_r*cylinder_r);
    
    auto center = glm::vec3(0.0);
    entities.reserve(entities.size() + data.num_atoms);
    for(int i = 0; i < data.num_atoms; ++i){
        auto &a = data.atoms[i];
        auto m_e = new MeshEntity(i);

        m_e->albedo = glm::normalize(symbol_to_color_lut[std::string(a.symbol)]);

        m_e->mesh = &graphics::sphere;
        m_e->scale = glm::mat3(sphere_r);
        m_e->position = a.position;

        entities.push_back(m_e);

        center += a.position;
    }
    center /= data.num_atoms;

    // Alternate method for reserving and also incase double bonds are written as two singles
    //int num_bond_entities;
    //std::unordered_map<int, MolBondType> atom_to_bonds;
    //for(int i = 0; i < data.num_bonds; ++i){
    //    auto &b = data.bonds[i];
    //    auto type = b.type;
    //    if(type == MolBondType::UNKNOWN) type = MolBondType::SINGLE;

    //    // @note that there is a maximum of 999 atoms so bitwise or is valid for >32 bit int
    //    int key = b.atom_1_index | (b.atom_2_index << 16);
    //    auto loc = atom_to_bonds.find(key);
    //    if(loc == atom_to_bonds.end()){
    //        atom_to_bonds[key] = type;
    //    } else {
    //        auto old_type = loc->second;
    //        // @note not guaranteed to be correct if mol file tries to overwrite bond type
    //        atom_to_bonds[key] = (MolBondType)((unsigned int)old_type + (unsigned int)type);
    //    }
    //    num_bond_entities += (int)type;
    //}

    //entities.reserve(entities.size() + num_bond_entities);
    //for(const auto &i : atom_to_bonds){
    //    // Get bits 32 --> 16
    //    int atom_2_index = i.first >> 16;
    //    // Get bits 16 --> 0 by removing 32 --> 16 since rest is empty
    //    int atom_1_index = i.first ^ (atom_2_index << 16);

    //    auto &atom_1 = data.atoms[atom_1_index];
    //    auto &atom_2 = data.atoms[atom_2_index];

    //    const MolBondType &type = i.second;


    for(int i = 0; i < data.num_bonds; ++i){
        auto &b = data.bonds[i];
        auto type = b.type;
        auto atom_1 = data.atoms[b.atom_1_index];
        auto atom_2 = data.atoms[b.atom_2_index];

        printf("Bond %d --> %d with type %u\n", b.atom_1_index, b.atom_2_index, type);
        switch(type){
            case MolBondType::SINGLE:
            {
                auto delta = atom_2.position - atom_1.position;
                auto distance = glm::length(delta); 

                auto m_e_1 = new MeshEntity(entities.size());
                auto m_e_2 = new MeshEntity(entities.size() + 1);
                m_e_1->mesh = &graphics::cylinder;
                m_e_2->mesh = &graphics::cylinder;

                m_e_1->albedo = glm::normalize(symbol_to_color_lut[std::string(atom_1.symbol)]);
                m_e_2->albedo = glm::normalize(symbol_to_color_lut[std::string(atom_2.symbol)]);

                m_e_1->position = atom_1.position;
                m_e_2->position = atom_2.position;

                scaleMat3(m_e_1->scale, glm::vec3(distance/2.0,  cylinder_r, cylinder_r));
                scaleMat3(m_e_2->scale, glm::vec3(distance/2.0,  cylinder_r, cylinder_r));

                m_e_1->rotation = quatAlignAxisToDirection(glm::vec3(1,0,0),  delta);
                m_e_2->rotation = quatAlignAxisToDirection(glm::vec3(1,0,0), -delta);

                entities.push_back(m_e_1);
                entities.push_back(m_e_2);

                auto m_e = new MeshEntity(entities.size());
                m_e->mesh = &graphics::cylinder;
                break;
            }
            case MolBondType::DOUBLE:
            {

                auto delta = atom_2.position - atom_1.position;
                auto distance = glm::length(delta); 
                auto v_offset = anyPerpendicular(delta);

                for(int j = 0; j < 2; ++j){
                    auto m_e_1 = new MeshEntity(entities.size());
                    auto m_e_2 = new MeshEntity(entities.size() + 1);
                    m_e_1->mesh = &graphics::cylinder;
                    m_e_2->mesh = &graphics::cylinder;

                    m_e_1->albedo = glm::normalize(symbol_to_color_lut[std::string(atom_1.symbol)]);
                    m_e_2->albedo = glm::normalize(symbol_to_color_lut[std::string(atom_2.symbol)]);

                    m_e_1->position = atom_1.position + bond_gap_r*v_offset;
                    m_e_2->position = atom_2.position + bond_gap_r*v_offset;

                    scaleMat3(m_e_1->scale, glm::vec3(distance/2.0, cylinder_r, cylinder_r));
                    scaleMat3(m_e_2->scale, glm::vec3(distance/2.0, cylinder_r, cylinder_r));

                    m_e_1->rotation = quatAlignAxisToDirection(glm::vec3(1,0,0),  delta);
                    m_e_2->rotation = quatAlignAxisToDirection(glm::vec3(1,0,0), -delta);

                    v_offset *= -1.0;

                    entities.push_back(m_e_1);
                    entities.push_back(m_e_2);
                }
                break;
            }
            case MolBondType::TRIPLE:
            {
                auto delta = atom_2.position - atom_1.position;
                auto distance = glm::length(delta); 
                auto v_offset = glm::normalize(anyPerpendicular(delta));

                for(int j = 0; j < 3; ++j){
                    auto m_e_1 = new MeshEntity(entities.size());
                    auto m_e_2 = new MeshEntity(entities.size() + 1);
                    m_e_1->mesh = &graphics::cylinder;
                    m_e_2->mesh = &graphics::cylinder;

                    m_e_1->albedo = glm::normalize(symbol_to_color_lut[std::string(atom_1.symbol)]);
                    m_e_2->albedo = glm::normalize(symbol_to_color_lut[std::string(atom_2.symbol)]);

                    m_e_1->position = atom_1.position + bond_gap_r*v_offset;
                    m_e_2->position = atom_2.position + bond_gap_r*v_offset;

                    scaleMat3(m_e_1->scale, glm::vec3(distance/2.0, cylinder_r, cylinder_r));
                    scaleMat3(m_e_2->scale, glm::vec3(distance/2.0, cylinder_r, cylinder_r));

                    m_e_1->rotation = quatAlignAxisToDirection(glm::vec3(1,0,0),  delta);
                    m_e_2->rotation = quatAlignAxisToDirection(glm::vec3(1,0,0), -delta);

                    v_offset = v_offset*glm::angleAxis((float)(2.0/3.0 * PI), delta);

                    entities.push_back(m_e_1);
                    entities.push_back(m_e_2);
                }
                break;
            }
            default:
                fprintf(stderr, "Unhandled bond type %u\n", type);
                break;
        }
    }
    return center;
}

#include <cstring>
#include <unordered_map>

#include <glm/glm.hpp>
#include "assets.hpp"
#include "glm/detail/func_geometric.hpp"
#include "glm/gtc/quaternion.hpp"

#include "loader.hpp"
#include "entities.hpp"
#include "graphics.hpp"
#include "utilities.hpp"

void loadMolFile(MolFile &data, std::string path){
    FILE *f;
    f=fopen(path.c_str(), "r");
    if (!f) {
        fprintf(stderr, "Error in reading mol file %s.\n", path.c_str());
        return;
    }

    printf("----------------Loading Mol File %s----------------\n", path.c_str());

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

static const float sphere_r = 0.6;
static const float cylinder_r = 0.035;
static const float bond_gap_r = 0.08;

glm::vec3 getColorFromSymbol(std::string symbol) {
    auto color_1_lu = symbol_to_color_lut.find(std::string(symbol));
    if(color_1_lu == symbol_to_color_lut.end()) return glm::normalize(glm::vec3(1.0));
    else                                        return glm::normalize(color_1_lu->second);
}

void createSingleBondEntities(glm::vec3 &pos_1, glm::vec3 &col_1, glm::vec3 &pos_2, glm::vec3 &col_2, std::vector<Entity*> &entities) {
    auto delta = pos_2 - pos_1;
    auto distance = glm::length(delta); 

    auto m_e_1 = new MeshEntity(entities.size());
    auto m_e_2 = new MeshEntity(entities.size() + 1);
    m_e_1->mesh = &graphics::cylinder;
    m_e_2->mesh = &graphics::cylinder;

    m_e_1->albedo = col_1;
    m_e_2->albedo = col_2;

    m_e_1->position = pos_1;
    m_e_2->position = pos_2;

    scaleMat3(m_e_1->scale, glm::vec3(distance/2.0,  cylinder_r, cylinder_r));
    scaleMat3(m_e_2->scale, glm::vec3(distance/2.0,  cylinder_r, cylinder_r));

    m_e_1->rotation = quatAlignAxisToDirection(glm::vec3(1,0,0),  delta);
    m_e_2->rotation = quatAlignAxisToDirection(glm::vec3(1,0,0), -delta);

    entities.push_back(m_e_1);
    entities.push_back(m_e_2);

    auto m_e = new MeshEntity(entities.size());
    m_e->mesh = &graphics::cylinder;
}

void createDoubleBondEntities(glm::vec3 &pos_1, glm::vec3 &col_1, glm::vec3 &pos_2, glm::vec3 &col_2, std::vector<Entity*> &entities) {
    auto v_offset = anyPerpendicular(pos_2 - pos_1);
    for(int j = 0; j < 2; ++j){
        auto offset_1 = pos_1 + bond_gap_r*v_offset;
        auto offset_2 = pos_2 + bond_gap_r*v_offset;
        createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities);
        v_offset *= -1.0;
    }
}

void createTripleBondEntities(glm::vec3 &pos_1, glm::vec3 &col_1, glm::vec3 &pos_2, glm::vec3 &col_2, std::vector<Entity*> &entities) {
    auto v_offset = anyPerpendicular(pos_2 - pos_1);
    for(int j = 0; j < 2; ++j){
        auto offset_1 = pos_1 + bond_gap_r*v_offset;
        auto offset_2 = pos_2 + bond_gap_r*v_offset;
        createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities);
        v_offset = v_offset*glm::angleAxis((float)(2.0/3.0 * PI), glm::normalize(pos_2 - pos_1));
    }
}

void createAtomEntity(glm::vec3 &pos, glm::vec3 &col, std::vector<Entity*> &entities) {
    auto m_e = new MeshEntity(entities.size());

    m_e->albedo = col;
    m_e->mesh = &graphics::sphere;
    m_e->scale = glm::mat3(sphere_r);
    m_e->position = pos;

    entities.push_back(m_e);
}

glm::vec3 createEntitiesFromMolFile(std::vector<Entity*> &entities, MolFile &data){
    auto center = glm::vec3(0.0);
    entities.reserve(entities.size() + data.num_atoms);
    for(int i = 0; i < data.num_atoms; ++i){
        auto &a = data.atoms[i];
        auto color = getColorFromSymbol(a.symbol);
        createAtomEntity(a.position, color, entities);
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

        auto color_1 = getColorFromSymbol(atom_1.symbol);
        auto color_2 = getColorFromSymbol(atom_2.symbol);

        printf("Bond %d --> %d with type %u\n", b.atom_1_index, b.atom_2_index, type);
        switch(type){
            case MolBondType::SINGLE:
            {
                createSingleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities);
                break;
            }
            case MolBondType::DOUBLE:
            {

                createDoubleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities);
                break;
            }
            case MolBondType::TRIPLE:
            {
                createTripleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities);
                break;
            }
            default:
                fprintf(stderr, "Unhandled bond type %u\n", type);
                break;
        }
    }
    return center;
}


void loadPdbFile(PdbFile &data, std::string path){
    FILE *f;
    f=fopen(path.c_str(), "r");
    if (!f) {
        fprintf(stderr, "Error in reading pdb file %s.\n", path.c_str());
        return;
    }

    printf("----------------Loading Pdb File %s----------------\n", path.c_str());

    // Assumes line is less than 1024 characters
    char line[1024];
    char record_name[7];

    data.polymer_models.emplace_back();
    const int default_polymer = data.polymer_models.size() - 1;
    int current_polymer = default_polymer;
    while(1) {
        if(fgets(line, 1024, f) == NULL) break;

        int len = strlen(line);
        if(len >= 6) {
            memcpy(record_name, &line[0], 6);
            record_name[6] = '\0';

            if(!strcmp(record_name, "MODEL")) {
                auto &model = data.polymer_models.emplace_back();
                sscanf(line, "MODEL %d", &model.serial);
                current_polymer = model.serial;
                printf("MODEL %d\n", current_polymer);
            } else if(!strcmp(record_name, "HETATM") || !strcmp(record_name, "ATOM  ")) {
                PdbAtom *atom;
                if(!strcmp(record_name, "HETATM") || current_polymer == -1) {
                    printf("HETATM: ");
                    atom = &data.heterogen_model.atoms.emplace_back();
                } else {
                    printf("ATOM  : ");
                    atom = &data.polymer_models[current_polymer].atoms.emplace_back();
                }
                sscanf(line, "HETATM%5d", &atom->serial);
                printf("Serial %d ", atom->serial);

                char pos_buf[25];
                memcpy(pos_buf, &line[30], 24);
                pos_buf[24] = '\0';
                sscanf(pos_buf, "%f %f %f", &atom->position.x, &atom->position.y, &atom->position.z);

                printf("Position %f %f %f ", atom->position.x, atom->position.y, atom->position.z);

                sscanf(&line[76], "%2s", atom->symbol);
                atom->symbol[1] = tolower(atom->symbol[1]);
                atom->symbol[2] = '\0';
                printf("Symbol %s\n", atom->symbol);
            } else if(!strcmp(record_name, "ENDMDL")) {
                current_polymer = default_polymer;
            }
        }
    }
}

glm::vec3 createEntitiesFromPdbFile(std::vector<Entity*> &entities, PdbFile &data){
    auto center = glm::vec3(0.0);

    entities.reserve(entities.size() + data.heterogen_model.atoms.size());
    for(auto &a : data.heterogen_model.atoms) {
        auto color = getColorFromSymbol(a.symbol);
        createAtomEntity(a.position, color, entities);
        center += a.position;
    }
    center /= data.heterogen_model.atoms.size();

    for(auto &m : data.polymer_models) {
        entities.reserve(entities.size() + m.atoms.size());
        for(auto &a : m.atoms) {
            auto color = getColorFromSymbol(a.symbol);
            createAtomEntity(a.position, color, entities);
            center += a.position;
        }
        center /= m.atoms.size();
    }

    return center;
}

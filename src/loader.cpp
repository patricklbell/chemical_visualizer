#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <set>

#include <glm/glm.hpp>
#include "controls.hpp"
#include "glm/detail/func_geometric.hpp"
#include "glm/gtc/quaternion.hpp"

#include "assets.hpp"
#include "loader.hpp"
#include "entities.hpp"
#include "graphics.hpp"
#include "utilities.hpp"

//
// General helper functions used by all loaders
//
static const float sphere_r = 0.7;
static const float cylinder_r = 0.2;
static const float bond_gap_r = 0.1;

glm::vec4 getColorFromSymbol(std::string symbol) {
    auto color_1_lu = symbol_to_color_lut.find(std::string(symbol));
    if (color_1_lu == symbol_to_color_lut.end())
    {
        fprintf(stderr, "Unknown element symbol %s", symbol.c_str());
        return glm::vec4(1.0);
    }
    else
    {
        return glm::vec4(color_1_lu->second / 255.f, 1.0);
    }
}

int createAtomEntity(const glm::vec3& pos, const glm::vec4& col, Entities& entities, float radius = sphere_r) {
    int id = entities.instanced_entities.size();
    auto& m_e = entities.instanced_entities.emplace_back();

    m_e.albedo = col;
    m_e.instance_mesh = InstancedMeshType::SPHERE;
    m_e.scale = glm::mat3(radius);
    m_e.position = pos;

    return id;
}

int createCylinderEntity(const glm::vec3& pos_1, const glm::vec3& pos_2, const glm::vec4& col, Entities& entities, float radius = cylinder_r) {
    auto delta = pos_2 - pos_1;
    auto distance = glm::length(delta);

    int id = entities.instanced_entities.size();
    auto& e = entities.instanced_entities.emplace_back();
    e.instance_mesh = InstancedMeshType::CYLINDER;
    scaleMat3(e.scale, glm::vec3(distance, radius, radius));

    e.position = pos_1;
    e.albedo = col;
    e.rotation = quatAlignAxisToDirection(glm::vec3(1, 0, 0), delta);

    return id;
}

int createSingleBondEntities(const glm::vec3& pos_1, const glm::vec4& col_1, const glm::vec3& pos_2, const glm::vec4& col_2, Entities& entities, float radius = cylinder_r) {
    auto delta = pos_2 - pos_1;
    auto distance = glm::length(delta);

    int id;
    for (int i = 0; i < 2; ++i) {
        id = entities.instanced_entities.size();
        auto& e = entities.instanced_entities.emplace_back();
        e.instance_mesh = InstancedMeshType::CYLINDER;
        scaleMat3(e.scale, glm::vec3(distance / 2.0, radius, radius));

        if (i == 0) {
            e.position = pos_1;
            e.albedo = col_1;
            e.rotation = quatAlignAxisToDirection(glm::vec3(1, 0, 0), delta);
        }
        else {
            e.position = pos_2;
            e.albedo = col_2;
            e.rotation = quatAlignAxisToDirection(glm::vec3(1, 0, 0), -delta);
        }
    }
    return id;
}

void createDebugCartesian(const glm::vec3& p, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, Entities& entities, float r = cylinder_r) {
    createAtomEntity(p, glm::vec4(1.0, 1.0, 1.0, 0.5), entities, r);
    createCylinderEntity(p, p + a, glm::vec4(1, 0, 0, 0.5), entities, r);
    createCylinderEntity(p, p + b, glm::vec4(0, 1, 0, 0.5), entities, r);
    createCylinderEntity(p, p + c, glm::vec4(0, 0, 1, 0.5), entities, r);
}

int createDoubleBondEntities(const glm::vec3& pos_1, const glm::vec4& col_1, const glm::vec3& pos_2, const glm::vec4& col_2, Entities& entities, float radius = cylinder_r) {
    auto v_offset = anyPerpendicular(pos_2 - pos_1);
    int id;
    for (int j = 0; j < 2; ++j) {
        auto offset_1 = pos_1 + bond_gap_r * v_offset;
        auto offset_2 = pos_2 + bond_gap_r * v_offset;
        id = createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities, radius / 2.0);
        v_offset *= -1.0;
    }
    return id;
}

int createTripleBondEntities(const glm::vec3& pos_1, const glm::vec4& col_1, const glm::vec3& pos_2, const glm::vec4& col_2, Entities& entities, float radius = cylinder_r) {
    auto v_offset = anyPerpendicular(pos_2 - pos_1);
    int id;
    for (int j = 0; j < 3; ++j) {
        auto offset_1 = pos_1 + bond_gap_r * v_offset;
        auto offset_2 = pos_2 + bond_gap_r * v_offset;
        id = createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities, radius / 2.0);
        v_offset = v_offset * glm::angleAxis((float)(2.0 / 3.0 * glm::pi<float>()), glm::normalize(pos_2 - pos_1));
    }
    return id;
}

//
// Mol file loading and entity creation
//

// @todo error checking
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

    int num_atoms, num_bonds;
    fscanf(f, "%d %d %*[^\n] ", &num_atoms, &num_bonds);
    printf("Num of atoms is %d and bonds is %d\n", num_atoms, num_bonds);

    data.atoms.resize(num_atoms);
    for(int i = 0; i < num_atoms; ++i){
        auto &a = data.atoms[i];

        fscanf(f, "%f %f %f", &a.position.x, &a.position.y, &a.position.z);
        fscanf(f, "%s", a.symbol);
        int md;
        fscanf(f, "%d %*d %d %*d %d %*[^\n]", &md, &a.hydrogen_count, &a.valence);
        a.mass_difference = (float)md;

        printf("Atom %s position is x: %f y: %f z: %f\n", a.symbol, a.position.x, a.position.y, a.position.z);
    }

    data.bonds.resize(num_bonds);
    for(int i = 0; i < num_bonds; ++i){
        auto &b = data.bonds[i];

        fscanf(f, "%d %d %u %*[^\n]", &b.atom_1_index, &b.atom_2_index, &b.type);
        b.atom_1_index--;
        b.atom_2_index--;
        printf("Bond %d from atom %d --> %d\n", i, b.atom_1_index, b.atom_2_index);
    }
}

void MolFile::clear() {
    atoms.clear();
    bonds.clear();
}

void createEntitiesFromMolFile(Entities &entities, MolFile &data, Camera &camera){
    static const float relative_cylinder_size = 0.1;
    static const float relative_sphere_size   = 0.26;

    float avg_bond_length = 0.0;
    for(auto &b : data.bonds){
        auto type = b.type;
        auto &atom_1 = data.atoms[b.atom_1_index];
        auto &atom_2 = data.atoms[b.atom_2_index];
        avg_bond_length += glm::length(atom_2.position - atom_1.position);
    }
    avg_bond_length /= data.bonds.size();
    float calc_cylinder_r = avg_bond_length*relative_cylinder_size;
    float calc_sphere_r   = avg_bond_length*relative_sphere_size;

    auto center = glm::vec3(0.0);
    entities.instanced_entities.reserve(entities.instanced_entities.size() + data.atoms.size());
    for(auto &a : data.atoms){
        auto color = getColorFromSymbol(a.symbol);
        createAtomEntity(a.position, color, entities, calc_sphere_r);
        center += a.position;
    }
    center /= data.atoms.size();

    float max_distance = 0.0;
    for (auto& a : data.atoms) {
        auto dist = glm::length(a.position - center);
        if (dist > max_distance) {
            max_distance = dist;
        }
    }

    camera.selected = true;
    camera.selection_position = center;
    camera.selection_radius = max_distance;

    // @alternate method for reserving and also incase double bonds are written as two singles
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
    
    // Each bond is at least a single bond which is two entities
    entities.instanced_entities.reserve(entities.instanced_entities.size() + data.bonds.size()*2);
    for(auto& b : data.bonds) {
        auto type = b.type;
        auto atom_1 = data.atoms[b.atom_1_index];
        auto atom_2 = data.atoms[b.atom_2_index];

        auto color_1 = getColorFromSymbol(atom_1.symbol);
        auto color_2 = getColorFromSymbol(atom_2.symbol);

        printf("Bond %d --> %d with type %u\n", b.atom_1_index, b.atom_2_index, type);
        switch(type){
            case MolBondType::SINGLE:
            {
                createSingleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            case MolBondType::DOUBLE:
            {

                createDoubleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            case MolBondType::TRIPLE:
            {
                createTripleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            default:
                fprintf(stderr, "Unhandled bond type %u\n", type);
                break;
        }
    }
}

//
// Pdb file loading and entity creation
//
long hashBondPair(int id_1, int id_2) {
    if(id_2 < id_1) return (long)id_1 ^ ((long)id_2 << 16);
    else            return (long)id_2 ^ ((long)id_1 << 16);
}

int hashResidueSeqChain(int res_seq, char chain_id) {
    // Since res <= 9999
    return res_seq + 100000*chain_id;
}

void createResidueId(int seq_num, char res_name[4], char chain_id, char i_code, char result[10]) {
    sprintf(result, "%3s%c%4d%c", res_name, chain_id, seq_num, i_code); 
}

glm::vec4 getColorFromResidueName(std::string residue) {
    auto color_1_lu = residue_to_color_lut.find(std::string(residue));
    if (color_1_lu == residue_to_color_lut.end()) return glm::vec4(unknown_residue_color, 1.0);
    else                                        return glm::vec4(color_1_lu->second, 1.0);
}

// Dictionary that provides extra information about residue bonds,
// Its format is not the same as a normal pdb file
void loadPdbDictionaryFile(PdbDictionary &dict, std::string_view path) {
    FILE *f;
    f=fopen(path.data(), "r");
    if (!f) {
        fprintf(stderr, "Error in reading pdb dictionary file %s.\n", path.data());
        return;
    }

    printf("----------------Loading Pdb Dictionary File %s----------------\n", path.data());

    // @note Assumes line is less than 1024 characters
    char line[1024];
    char record_name[8];

    std::unordered_map<std::string, PdbDictionaryConnect> *current_residue_connections = nullptr;
    while(1) {
        if(fgets(line, 1024, f) == NULL) break;

        if(strlen(line) >= 6) {
            substrString(line, 0, 6, record_name);

            if(!strcmp(record_name, "RESIDUE")) {
                char res_name[4];
                substrString(line, 10, 12,  res_name);

                current_residue_connections = &dict.residues.try_emplace(std::string(res_name)).first->second;

                //printf("RESDUE: res_name %s \n", res_name); // @debug
            } else if(!strcmp(record_name, "CONECT ")) { 
                if(current_residue_connections == nullptr) continue;

                char atom_name[5];
                substrString(line, 12, 15, atom_name);

                int num_bonds = 0;
                substrInt(line, 19, 19, &num_bonds);

                char bonded_atom_name[5];

                //printf("CONECT: %s num %d --> ", atom_name, num_bonds); // @debug

                for(int i = 0; i < num_bonds; i++) {
                    substrString(line, 21 + i*5, 24 + i*5, bonded_atom_name);

                    //if(i == num_bonds - 1)
                    //    printf("%s\n", bonded_atom_name);
                    //else
                    //    printf("%s, ", bonded_atom_name);

                    // @todo make this hashing system better
                    std::string atom_name_1 = std::string(atom_name);
                    std::string atom_name_2 = std::string(bonded_atom_name);
                    std::string name_hash;
                    if(atom_name_1.compare(atom_name_2) == -1) // Compares alphabetically
                        name_hash = atom_name_1 + atom_name_2;
                    else
                        name_hash = atom_name_2 + atom_name_1;
                    auto lu = current_residue_connections->find(name_hash);

                    // Check if a connection exists, if not, create one
                    PdbDictionaryConnect *connect;
                    if(lu == current_residue_connections->end()) {
                        connect = &current_residue_connections->try_emplace(name_hash).first->second;
                        connect->atom_name_1 = atom_name_1;
                        connect->atom_name_2 = atom_name_2;
                        connect->type = PdbConnectionType::UNKNOWN; // 0
                    } else {
                        connect = &lu->second;
                    }
                    connect->type = (PdbConnectionType)glm::min((unsigned int)connect->type + 1, (unsigned int)PdbConnectionType::TRIPLE_REDUNDANT); // Fix this C++!
                }
            }
        }
    }
}

void loadPdbFile(PdbFile &data, std::string path, PdbDictionary *dict){
    FILE *f;
    f=fopen(path.c_str(), "r");
    if (!f) {
        fprintf(stderr, "Error in reading pdb file %s.\n", path.c_str());
        return;
    }

    printf("----------------Loading Pdb File %s----------------\n", path.c_str());

    // @note Assumes line is less than 1024 characters
    char line[1024];
    char record_name[7];

    auto model = &data.models.try_emplace(0).first->second;
    model->serial = 0;
    // Handle model command being missing
    bool first_model_flag = true;

    while(1) {
        if(fgets(line, 1024, f) == NULL) break;

        if(strlen(line) >= 6) {
            substrString(line, 0, 5, record_name);

            if(!strcmp(record_name, "MODEL ")) {
                int serial;
                sscanf(line, "MODEL %d", &serial);

                if(first_model_flag) {
                    first_model_flag = false;
                } else {
                    if (data.models.find(serial) != data.models.end()) {
                        fprintf(stderr, "Duplicate model serial number, %d\n", serial);
                        continue;
                    }
                    model = &data.models.try_emplace(serial).first->second;
                }
                
                model->serial = serial;
                printf("MODEL %d\n", model->serial);
            } else if(!strcmp(record_name, "HETATM") || !strcmp(record_name, "ATOM  ")) {
                int serial;
                substrSscanf(line, 6, 10, " %d", &serial);

                // @note this overwrites duplicate ids
                auto lu0 = model->atoms.find(serial);
                if(lu0 != model->atoms.end()) printf("Collision for serial %d\n", serial);
                
                PdbAtom &atom = model->atoms[serial];
                atom.serial = serial;
                atom.is_heterogen = !strcmp(record_name, "HETATM");

                substrString(line, 13, 16,  atom.name);

                substrString(line, 17, 19,  atom.res_name);
                substrChar  (line, 21,     &atom.chain_id);
                substrInt   (line, 22, 25, &atom.res_seq);

                auto residue_id = hashResidueSeqChain(atom.res_seq, atom.chain_id);
                auto res_lu = model->residues.find(residue_id);
                // If residue doesn't exist create it
                if(res_lu == model->residues.end()) {
                    auto &residue = model->residues.try_emplace(residue_id).first->second;

                    memcpy(residue.res_name, atom.res_name, 4);
                    residue.chain_id = atom.chain_id;
                    residue.i_code   = atom.i_code;
                    residue.res_seq  = atom.res_seq;

                    // Add new residue to chain
                    auto chain_lu = model->chains.find(atom.chain_id);
                    // If chain doesn't exist create it
                    // @note doesn't account for empty chain_id i.e ' '
                    if(chain_lu == model->chains.end()) {
                        auto &chain = model->chains.try_emplace(atom.chain_id).first->second;
                        chain.chain_id = atom.chain_id;
                        // @note pointers to values in unordered_map are stable 
                        chain.residues.push_back(&residue);

                        printf("CHAIN : %c\n", chain.chain_id); // @debug
                    } else {
                        auto &chain = chain_lu->second;
                        chain.residues.push_back(&residue);
                    }

                    residue.atom_ids.push_back(atom.serial);
                    residue.atom_name_id[atom.name] = atom.serial;

                    printf("RESDUE: res_name %s chain_id %c i_code %c res_seq %d\n", 
                            residue.res_name, residue.chain_id, residue.i_code, residue.res_seq); // @debug
                } else {
                    auto &residue = res_lu->second;
                    residue.atom_ids.push_back(atom.serial);
                    residue.atom_name_id[atom.name] = atom.serial;
                }

                substrFloat(line, 30, 37, &atom.position.x);
                substrFloat(line, 38, 45, &atom.position.y);
                substrFloat(line, 46, 53, &atom.position.z);
                sscanf(&line[76], "%2s", atom.symbol);
                atom.symbol[2] = '\0';

                //if(atom.is_heterogen) printf("HETATM:");
                //else                  printf("ATOM  :");
                //printf(" serial %4d name %s symbol %s res_name %s chain_id %c i_code %c res_seq %3d\n", 
                //        atom.serial, atom.name, atom.symbol, atom.res_name, atom.chain_id, atom.i_code, atom.res_seq);
            } else if(!strcmp(record_name, "CONECT")) { 
                int atom_id;
                int connect_ids[4];

                substrInt(line, 6, 10, &atom_id);
                int matches = 0;
                for(; matches < 3; matches++) {
                    if(substrSscanf(line, 11 + 5*matches, 11 + 5*matches + 4, " %d", &connect_ids[matches]) != 1) 
                        break; 
                }

                printf("CONECT: %d --> ", atom_id);
                for(int i = 0; i < matches; ++i) {
                    printf("%d, ", connect_ids[i]);
                    auto hash = hashBondPair(atom_id, connect_ids[i]);
                    auto lu = model->connections.find(hash);
                    if(lu != model->connections.end()) {
                        auto &b = lu->second;
                        b.type = PdbConnectionType((unsigned int)b.type + 1);
                    } else {
                        auto b = &model->connections[hash];
                        b->atom_1_id = atom_id;
                        b->atom_2_id = connect_ids[i];
                        b->type = PdbConnectionType::SINGLE;
                    }
                }
                printf("\n");
            } else if(!strcmp(record_name, "HELIX ")) { 
                int serial;
                substrInt(line, 7, 9, &serial);
                auto &helix = model->helices.try_emplace(serial).first->second;
                helix.serial = serial;
                
                substrString(line, 11, 13,  helix.id_code);

                substrString(line, 15, 17,  helix.init_res_name);
                substrChar  (line, 19,     &helix.init_chain_id);
                substrInt   (line, 21, 24, &helix.init_seq_num);
                substrChar  (line, 25,     &helix.init_i_code);
                
                substrString(line, 27, 29,  helix.end_res_name);
                substrChar  (line, 31,     &helix.end_chain_id);
                substrInt   (line, 33, 36, &helix.end_seq_num);
                substrChar  (line, 37,     &helix.end_i_code);

                substrInt   (line, 71, 75, &helix.seq_length);

                printf("HELIX : id_code %s init_seq_num %d end_seq_num %d seq_length %d\n",
                        helix.id_code, helix.init_seq_num, helix.end_seq_num, helix.seq_length); // @debug
            } else if(!strcmp(record_name, "SHEET ")) { 
                char sheet_id[4];
                substrString(line, 11, 13, sheet_id);

                bool is_first = false;
                auto lu = model->sheets.find(sheet_id);
                PdbSheet *sheet;
                if(lu == model->sheets.end()) {
                    // Create and modify sheet only first time
                    sheet = &model->sheets.try_emplace(sheet_id).first->second;
                    memcpy(sheet->sheet_id, sheet_id, 4);
                    substrInt(line, 14, 15, &sheet->num_strands);
                    is_first = true;
                    printf("SHEET : sheet_id %s num_strands %d\n",
                            sheet->sheet_id, sheet->num_strands); // @debug
                } else {
                    sheet = &model->sheets[sheet_id];
                }
                int strand_id;
                substrInt  (line,  7,  9,  &strand_id);
                auto &strand = sheet->strands.try_emplace(strand_id).first->second;
                strand.strand_id = strand_id;
                strand.is_first  = is_first;
                
                substrString(line, 17, 19,  strand.init_res_name);
                substrChar  (line, 21,     &strand.init_chain_id);
                substrInt   (line, 22, 25, &strand.init_seq_num);
                substrChar  (line, 26,     &strand.init_i_code);
                
                substrString(line, 28, 30,  strand.end_res_name);
                substrChar  (line, 32,     &strand.end_chain_id);
                substrInt   (line, 33, 36, &strand.end_seq_num);
                substrChar  (line, 37,     &strand.end_i_code);

                if(!is_first) {
                    substrString(line, 41, 44,  strand.cur_atom_name);
                    substrString(line, 45, 47,  strand.cur_res_name);
                    substrChar  (line, 49,     &strand.cur_chain_id);
                    substrInt   (line, 50, 53, &strand.cur_seq_num);
                    substrChar  (line, 54,     &strand.cur_i_code);
                    
                    substrString(line, 56, 59,  strand.prev_atom_name);
                    substrString(line, 60, 62,  strand.prev_res_name);
                    substrChar  (line, 64,     &strand.prev_chain_id);
                    substrInt   (line, 65, 68, &strand.prev_seq_num);
                    substrChar  (line, 69,     &strand.prev_i_code);

                } 
                printf("STRAND: sheet_id %s init_res_name %s init_seq_num %3d end_res_name %s end_seq_num %3d\n",
                        sheet->sheet_id, strand.init_res_name, strand.init_seq_num, strand.end_res_name, strand.end_seq_num); // @debug
            } else if(!strcmp(record_name, "ENDMDL")) {
                printf("ENDMDL\n");
            }
        }
    }

    // Processing after loading from file
    for(auto &pair : data.models) {
        auto& m = pair.second;
        // Label residue type based on secondary structures
        for(auto &j : m.helices) {
            auto &helix = j.second;

            if(helix.init_chain_id != helix.end_chain_id) fprintf(stderr, "Error: helix chain init and end not equal");
            char chain_id = helix.init_chain_id;

            for(int seq_num = helix.init_seq_num; seq_num <= helix.end_seq_num; seq_num++) {
                auto residue_id = hashResidueSeqChain(seq_num, chain_id);
                auto lu = m.residues.find(residue_id);
                if(lu == m.residues.end()) continue;

                auto &residue = lu->second;
                residue.type = PdbResidueType::HELIX;
            }
        }
        for(auto &j : m.sheets) {
            auto &sheet = j.second;
            
            for(auto &k : sheet.strands) {
                auto &strand = k.second;

                if(strand.init_chain_id != strand.end_chain_id) fprintf(stderr, "Error: strand chain_id init and end not equal");
                char chain_id = strand.init_chain_id;

                for(int seq_num = strand.init_seq_num; seq_num <= strand.end_seq_num; seq_num++) {
                    auto residue_id = hashResidueSeqChain(seq_num, chain_id);
                    auto lu = m.residues.find(residue_id);
                    if(lu == m.residues.end()) continue;

                    auto &residue = lu->second;
                    residue.type = PdbResidueType::STRAND;
                }
            }
        }

        // Create connections for unspecified residues from dictionary
        if(dict != nullptr) {
            for(auto &rp : m.residues) {
                auto &residue = rp.second;

                auto res_lu = dict->residues.find(std::string(residue.res_name));

                // The residue isn't in our dictionary
                if(res_lu == dict->residues.end()) {
                    fprintf(stderr, "RESIDUE: res_name %s chain_id %c i_code %c res_seq %d NOT IN DICTIONARY, Skipping\n", 
                            residue.res_name, residue.chain_id, residue.i_code, residue.res_seq); // @debug
                    continue;
                }

                auto &residue_connections = res_lu->second;
                for(auto &ap : residue_connections) {
                    auto atom_1_lu = residue.atom_name_id.find(ap.second.atom_name_1);
                    if(atom_1_lu == residue.atom_name_id.end()) {
                        //fprintf(stderr, "Atom %s in residue %s's dictionary doesn't exist in file\n", ap.second.atom_name_1.c_str(), residue.res_name);
                        continue;
                    }
                    auto &atom_1_id = atom_1_lu->second;
                    auto atom_2_lu = residue.atom_name_id.find(ap.second.atom_name_2);
                    if(atom_2_lu == residue.atom_name_id.end()) {
                        //fprintf(stderr, "Atom %s in residue %s's dictionary doesn't exist in file\n", ap.second.atom_name_2.c_str(), residue.res_name);
                        continue;
                    }
                    auto &atom_2_id = atom_2_lu->second;

                    auto hash = hashBondPair(atom_1_id, atom_2_id);
                    auto lu = m.connections.find(hash);
                    if(lu == m.connections.end()) {
                        //printf("DICT CONECT: %d --> %d of type %d\n", atom_1_id, atom_2_id, ap.second.type);
                        auto &b = m.connections.try_emplace(hash).first->second;
                        b.atom_1_id = atom_1_id;
                        b.atom_2_id = atom_2_id;
                        b.type = ap.second.type;
                    }
                }
            }
        }
    }
}

// Mesh which blends between different secondary structures in one polypeptide chain 
void createPolypeptideChainsMeshes(PdbChain &chain, std::vector<PeptidePlane> &planes) {
    constexpr float r = 0.25;
    // By default colors based on secondary structures
    constexpr bool do_chain_colors = true;
    constexpr bool do_residue_colors = false;

    // @todo these should probably be dynamic, @note num_splines has to be multiples of 8
    constexpr int num_splines = 16, num_points_per_spline = 12;
    glm::vec2 circle_profile[num_splines];
    glm::vec2 circle_profile_normals[num_splines];
    createCircularProfile(num_splines, circle_profile, r);
    createCircularProfileNormals(num_splines, circle_profile_normals);

    glm::vec2 ellipse_profile[num_splines];
    glm::vec2 ellipse_profile_normals[num_splines];
    createEllipseProfile(num_splines, ellipse_profile, 3.f*r, r/2.f);
    createEllipseProfileNormals(num_splines, ellipse_profile_normals, 3.f*r, r/2.f);

    glm::vec2 rectangle_profile[num_splines];
    glm::vec2 rectangle_profile_normals[num_splines];
    createRectangleProfile(num_splines, rectangle_profile, 6.f*r, 1.5f*r);
    createRectangleProfileNormals(num_splines, rectangle_profile_normals, 6.f*r, 1.5f*r);

    glm::vec2 ribbon_profile[num_splines];
    glm::vec2 ribbon_profile_normals[num_splines];
    createRectangleProfile(num_splines, ribbon_profile, 5.f*r, 0.5f*r);
    createRectangleProfileNormals(num_splines, ribbon_profile_normals, 5.f*r, 0.5f*r);

    // 4 point versions of rectangle profiles, used when not transitioning
    glm::vec2 rectangle_profile_simple[4];
    glm::vec2 rectangle_profile_simple_normals_doubled[8];
    createRectangleProfile(4, rectangle_profile_simple, 6.f * r, 1.5f * r);
    createRectangleProfileNormalsCorrect(rectangle_profile_simple_normals_doubled, 6.f * r, 1.5f * r);

    glm::vec2 ribbon_profile_simple[4];
    glm::vec2 ribbon_profile_simple_normals_doubled[8];
    createRectangleProfile(4, ribbon_profile_simple, 5.f * r, 0.5f * r);
    createRectangleProfileNormalsCorrect(ribbon_profile_simple_normals_doubled, 5.f * r, 0.5f * r);

    // [0,1] point to exavluate to get arrow head beginning 
    constexpr float arrow_length = 0.7;
    glm::vec2 arrow_profile [4] ;
    glm::vec2 arrow_profile_normals_doubled[8];
    createRectangleProfile(4, arrow_profile, 10.f * r, 1.5f * r);
    createRectangleProfileNormalsCorrect(arrow_profile_normals_doubled, 10.f * r, 1.5f * r);
    glm::vec2 arrow_tip_profile[4];
    glm::vec2 arrow_tip_profile_normals_doubled[8];
    createRectangleProfile(4, arrow_tip_profile, 0.f * r, 1.f * r);
    createRectangleProfileNormalsCorrect(arrow_tip_profile_normals_doubled, 0.f * r, 1.f * r);
    glm::vec2 arrow_transition_profile[num_splines];
    glm::vec2 arrow_transition_profile_normals[num_splines];
    createRectangleProfile(num_splines, arrow_transition_profile, 10.f * r, 1.5f * r);
    createRectangleProfileNormals(num_splines, arrow_transition_profile_normals, 10.f * r, 1.5f * r);

    // @todo Handle small number of planes
    if(planes.size() < 3) return;
    
    // Points describing spline tube's surface
    glm::vec3 spline_tube [num_points_per_spline][num_splines];
    glm::vec3 normals_tube[num_points_per_spline][num_splines];

    const int num_planes = planes.size();

    glm::vec3 previous_normal = planes[0].normal;

    for(int i = 0; i < num_planes - 1; i++) {
        auto& mesh = planes[i].residue_1->mesh;

        // @note middle peptide or similar might reduce floating point errors
        // need to offset all vertex positions then

        mesh.draw_mode = GL_TRIANGLES;
        mesh.draw_type = GL_UNSIGNED_SHORT;

        // @todo make parts of mesh different colors
        mesh.num_materials = 1;
        mesh.draw_start = (GLint*)malloc(sizeof(GLint) * mesh.num_materials);
        mesh.draw_count = (GLint*)malloc(sizeof(GLint) * mesh.num_materials);
        mesh.draw_start[0] = 0;
        mesh.draw_count[0] = 0;

        //createDebugCartesian(planes[i].position, 0.5f*planes[i].normal, 0.5f*planes[i].right, 0.5f*planes[i].forward, entities, 0.02);

        glm::vec2 *pf1, *pf2, *pfn1, *pfn2;
        // Flag set when rectangle profiles which aren't transitioning
        const bool do_simple_profile = (planes[i].residue_1->type == PdbResidueType::HELIX && planes[i + 1].residue_1->type == PdbResidueType::HELIX) ||
                                       (planes[i].residue_1->type == PdbResidueType::STRAND && planes[i + 1].residue_1->type == PdbResidueType::STRAND); 
        switch (planes[i].residue_1->type) {
            case PdbResidueType::COIL:
                {
                    pf1  = circle_profile;
                    pfn1 = circle_profile_normals;
                    break;
                }
            case PdbResidueType::HELIX:
                {
                    pf1  = ribbon_profile;
                    pfn1 = ribbon_profile_normals;
                    break;
                }
            case PdbResidueType::STRAND:
                {
                    pf1  = rectangle_profile;
                    pfn1 = rectangle_profile_normals;
                    break;
                }
            default:
                {
                    fprintf(stderr, "unknown residue type 1\n");
                    pf1  = circle_profile;
                    pfn1 = circle_profile_normals;
                    break;
                }
        }
        switch (planes[i+1].residue_1->type) {
            case PdbResidueType::COIL:
                {
                    pf2  = circle_profile;
                    pfn2 = circle_profile_normals;
                    break;
                }
            case PdbResidueType::HELIX:
                {
                    pf2  = ribbon_profile;
                    pfn2 = ribbon_profile_normals;
                    break;
                }
            case PdbResidueType::STRAND:
                {
                    pf2  = rectangle_profile;
                    pfn2 = rectangle_profile_normals;
                    break;
                }
            default:
                {
                    fprintf(stderr, "unknown residue type 2\n");
                    pf2  = circle_profile;
                    pfn2 = circle_profile_normals;
                    break;
                }
        }

        int p1_i = i-1;
        int p4_i = i+2;
        if(i == 0) {
            p1_i = i;
        } else if(i == num_planes - 2) {
            p4_i = i+1;
        }

        if(planes[i].residue_1->type == PdbResidueType::STRAND && planes[i+1].residue_1->type == PdbResidueType::COIL) {
            auto arrow_base_plane = createPartialHermiteSplineNormalsBetweenProfiles(num_points_per_spline, 4, 
                    rectangle_profile_simple, rectangle_profile_simple_normals_doubled, rectangle_profile_simple, rectangle_profile_simple_normals_doubled,
                    planes[p1_i], planes[i], planes[i+1], planes[p4_i], 
                    &spline_tube[0][0], &normals_tube[0][0], 
                    previous_normal, arrow_length, true);
            createClosedSurfaceFromSplinesDoubledNormals(4, num_points_per_spline, &spline_tube[0][0], &normals_tube[0][0], mesh);

            glm::vec3 projected_profile[4];
            projectPointsOnPlane(4, arrow_base_plane.position, arrow_base_plane.normal, arrow_base_plane.right, 
                    arrow_profile, projected_profile);
            createClosedFacedFromProfile(arrow_base_plane.position, 4, projected_profile, mesh, true);

            if(i == num_planes - 2) {
                // hack to make longer arrows at end points
                planes[i+1].position += arrow_base_plane.forward;
                createHermiteSplineNormalsBetweenProfiles(2, 4, 
                        arrow_profile, arrow_profile_normals_doubled, arrow_tip_profile, arrow_tip_profile_normals_doubled, 
                        planes[i], arrow_base_plane, planes[i+1], planes[p4_i], 
                        &spline_tube[0][0], &normals_tube[0][0], 
                        previous_normal, true);
                createClosedSurfaceFromSplinesDoubledNormals(4, 2, &spline_tube[0][0], &normals_tube[0][0], mesh);
            } else {
                createHermiteSplineNormalsBetweenProfiles(2, num_splines, 
                        arrow_transition_profile, arrow_transition_profile_normals, circle_profile, circle_profile_normals, 
                        planes[i], arrow_base_plane, planes[i+1], planes[p4_i],
                        &spline_tube[0][0], &normals_tube[0][0], 
                    previous_normal);
                createClosedSurfaceFromSplinesNormals(num_splines, 2, &spline_tube[0][0], &normals_tube[0][0], mesh);
            }

            //projectPointsOnPlane(num_splines, planes[i+1].position, planes[i+1].right, planes[i+1].normal, 
            //        circle_profile, projected_profile);
            //createClosedFacedFromProfile(planes[i+1].position, num_splines, projected_profile, mesh, true);
        } else {
            if (i == 0) {
                auto bref = (planes[i].right + planes[i+1].right) / 2.f;
                auto tn = glm::normalize(planes[i].forward);
                glm::vec3 bn;
                if (glm::length(previous_normal) > 0.001) {
                    bn = glm::normalize(glm::cross(tn, previous_normal));
                }
                else {
                    bn = perpendicularComponent(bref, tn);
                }
                auto n = glm::cross(bn, tn);

                glm::vec3 projected_profile[num_splines];
                projectPointsOnPlane(num_splines, planes[i].position, bn, n, pf1, projected_profile);
                createClosedFacedFromProfile(planes[i].position, num_splines, projected_profile, mesh, true);
            }

            // Rectange normals are wrong if they are smooth shaded so we need a special case
            if (do_simple_profile) {
                if (planes[i].residue_1->type == PdbResidueType::HELIX) {
                    createHermiteSplineNormalsBetweenProfiles(num_points_per_spline, 4,
                        ribbon_profile_simple, ribbon_profile_simple_normals_doubled, ribbon_profile_simple, ribbon_profile_simple_normals_doubled,
                        planes[p1_i], planes[i], planes[i + 1], planes[p4_i], 
                        &spline_tube[0][0], &normals_tube[0][0], 
                        previous_normal, true);
                }
                else if (planes[i].residue_1->type == PdbResidueType::STRAND) {
                    createHermiteSplineNormalsBetweenProfiles(num_points_per_spline, 4,
                        rectangle_profile_simple, rectangle_profile_simple_normals_doubled, rectangle_profile_simple, rectangle_profile_simple_normals_doubled,
                        planes[p1_i], planes[i], planes[i + 1], planes[p4_i], 
                        &spline_tube[0][0], &normals_tube[0][0], 
                        previous_normal, true);
                }
                createClosedSurfaceFromSplinesDoubledNormals(4, num_points_per_spline, &spline_tube[0][0], &normals_tube[0][0], mesh);
            }
            else {
                createHermiteSplineNormalsBetweenProfiles(num_points_per_spline, num_splines, pf1, pfn1, pf2, pfn2, 
                        planes[p1_i], planes[i], planes[i+1], planes[p4_i], &spline_tube[0][0], 
                        &normals_tube[0][0], previous_normal);
                createClosedSurfaceFromSplinesNormals(num_splines, num_points_per_spline, &spline_tube[0][0], &normals_tube[0][0], mesh);
            }

            if(i == num_planes - 2) {
                auto tn = glm::normalize(planes[i+1].forward);
                auto n = previous_normal;
                auto bn = glm::cross(tn, n);

                glm::vec3 projected_profile[num_splines];
                projectPointsOnPlane(num_splines, planes[i+1].position, bn, n, pf2, projected_profile);
                createClosedFacedFromProfile(planes[i+1].position, num_splines, projected_profile, mesh, false);
            }
        }
        createMeshVao(mesh);
        planes[i].residue_1->mesh_generated = true;
    }
}

void createPdbModelMeshes(PdbModel &model) {
    for (auto& p : model.chains) {
        auto& chain = p.second;
        if (chain.residues.size() < 2) continue;

        std::vector<PeptidePlane> peptide_planes;

        //printf("Chain %c\n", chain.chain_id);
        int plane_i = 0;
        for (int i = 0; i < chain.residues.size() - 1; i++) {
            // Determine plane of connection between residues through peptide bond
            auto& r1 = chain.residues[i];
            auto& r2 = chain.residues[i + 1];

            // @note Peptide bond forms between carboxylic group and amino group.
            // PDB should guarantee that residues' label the alpha carbon as "CA"
            // we just use the "O" of the amine to determine the perpendicular plane 
            // between adjacent amino acids of protein.
            auto lu_CA1 = r1->atom_name_id.find("CA  ");
            if (lu_CA1 == r1->atom_name_id.end()) continue;
            auto lu_CA2 = r2->atom_name_id.find("CA  ");
            if (lu_CA2 == r2->atom_name_id.end()) continue;
            //auto lu_O1  = r1->atom_name_id.find(" O  ");
            //if(lu_O1   == r1->atom_name_id.end()) continue;

            auto CA1 = model.atoms[lu_CA1->second].position;
            auto CA2 = model.atoms[lu_CA2->second].position;
            //auto O1  = model.atoms[lu_O1 ->second].position;

            auto& plane = peptide_planes.emplace_back();
            plane.residue_1 = r1;
            plane.residue_2 = r2;

            plane.position = (CA1 + CA2) / 2.f;
            plane.forward = glm::normalize(CA2 - CA1);
            plane.normal = glm::vec3(0);

            // This a more chemically accurate way of calulating plane frames but is less smooth
            //plane.forward = glm::normalize(CA2 - CA1);
            //plane.normal  = glm::normalize(glm::cross(plane.forward, O1 - CA1));
            //// To ensure perpendicularity calculate right
            //plane.right   = glm::normalize(glm::cross(plane.normal, plane.forward));

            //printf("Peptide Plane %d --> %d\n", i, i+1);
        }
        if (peptide_planes.size() == 0) continue;

        // Extends last peptide plane to pass through final alpha carbon,
        // @note I think this messes up the last frame so the closed face has incorrect normals
        /*{
            auto plane = peptide_planes[peptide_planes.size() - 1];
            auto lu_CA2 = plane.residue_2->atom_name_id.find("CA  ");
            if (lu_CA2 != plane.residue_2->atom_name_id.end()) {
                plane.position = model.atoms[lu_CA2->second].position;
                peptide_planes.emplace_back(plane);
            }
        }*/

        // https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Representations
        // Use Cardinal spline calculation of forward rather than using true plane for smoother result
        constexpr float c = 0.25;
        for (uint64_t i = 1; i < (uint64_t)peptide_planes.size() - 1; i++) {
            auto& p = peptide_planes;
            p[i].forward = (1 - c) * (p[i + 1].position - p[i - 1].position);
        }

        // Uses hermite spline derivatives to determine the normal so adjacent spline's orientations align
        for (uint64_t i = 0; i < (uint64_t)peptide_planes.size() - 1; i++) {
            auto& p = peptide_planes;
            auto n0 = -6.f * p[i].position - 4.f * p[i].forward + 6.f * p[i + 1].position - 2.f * p[i + 1].forward;
            n0 = perpendicularComponent(n0, p[i].forward);
            p[i].normal += n0 / 2.f;
            auto n1 = 6.f * p[i].position + 2.f * p[i].forward - 6.f * p[i + 1].position + 4.f * p[i + 1].forward;
            n1 = perpendicularComponent(n1, p[i + 1].forward);
            p[i + 1].normal += n1 / 2.f;
        }

        // Calulcate binormal and ensures adjacent normals aren't flipped
        for (uint64_t i = 0; i < peptide_planes.size(); i++) {
            auto& p = peptide_planes;
            p[i].normal = glm::normalize(p[i].normal);
            p[i].right = glm::normalize(glm::cross(p[i].forward, p[i].normal));
        }
        // This ensures the first and last planes normals are normalised since they only get half a normal
        {
            auto& pb = peptide_planes[0];
            auto& pe = peptide_planes[peptide_planes.size() - 1];

            pb.normal = glm::normalize(glm::cross(pb.forward, pb.right));
            pe.normal = glm::normalize(glm::cross(pe.forward, pe.right));
        }

        // Flip peptides that are >180 from previous
        for (uint64_t i = 1; i < peptide_planes.size(); i++) {
            if (glm::dot(peptide_planes[i].normal, peptide_planes[i - 1].normal) < 0) {
                peptide_planes[i].normal *= -1;
                peptide_planes[i].right *= -1;
            }
        }

        // Debug peptide planes
        /*for(int i = 1; i < peptide_planes.size(); i++) {
            createAtomEntity(peptide_planes[i].position, glm::vec4(1.0, 0.0, 0.0, 1.0), entities, calc_sphere_r);
            auto fwd = glm::normalize(peptide_planes[i].forward)*1.5f;
            createSingleBondEntities(peptide_planes[i].position - fwd*0.5f, glm::vec4(1,0,0,1), peptide_planes[i].position+fwd*0.5f, glm::vec4(1,0,0,1), entities,calc_cylinder_r/2.0f);
        }*/

        createPolypeptideChainsMeshes(chain, peptide_planes);
    }
}

// Way of updating entities in-place, definitely needs work @todo make this cleaner
void updateEntityColorsFromPdbModel(Entities& entities, PdbModel& model, PdbDrawSettings& settings) {
    for (const auto& p : model.connections) {
        auto& b = p.second;
        auto& type = b.type;

        auto atom_1_lu = model.atoms.find(b.atom_1_id);
        if (atom_1_lu == model.atoms.end()) continue;
        auto atom_1 = atom_1_lu->second;
        auto atom_2_lu = model.atoms.find(b.atom_2_id);
        if (atom_2_lu == model.atoms.end()) continue;
        auto atom_2 = atom_2_lu->second;

        auto color_1 = getColorFromSymbol(atom_1.symbol);
        auto color_2 = getColorFromSymbol(atom_2.symbol);

        bool is_heterogen_bond = atom_1.is_heterogen && atom_2.is_heterogen;

        if (is_heterogen_bond) {
            if (!settings.draw_hetero_atoms) continue;
            color_1.w = settings.hetero_atoms_alpha;
            color_2.w = settings.hetero_atoms_alpha;
        }
        else {
            if (!settings.draw_residue_atoms) continue;
            color_1.w = settings.residue_atoms_alpha;
            color_2.w = settings.residue_atoms_alpha;
        }

        if (b.entity_id != -1) {
            int num_entities = 0;
            switch (type) {
                case PdbConnectionType::SINGLE:
                case PdbConnectionType::SINGLE_REDUNDANT:
                {
                    num_entities = 2;
                    break;
                }
                case PdbConnectionType::DOUBLE:
                case PdbConnectionType::DOUBLE_REDUNDANT:
                {
                    num_entities = 4;
                    break;
                }
                case PdbConnectionType::TRIPLE:
                case PdbConnectionType::TRIPLE_REDUNDANT:
                {
                    num_entities = 6;
                    break;
                }
            }

            bool is_color1 = false;
            for (int id = b.entity_id; id > b.entity_id - num_entities && id >= 0 && id < entities.instanced_entities.size(); --id) {
                if (is_color1) entities.instanced_entities[id].albedo = color_1;
                else           entities.instanced_entities[id].albedo = color_2;
                is_color1 = !is_color1;
            }
        }
    }

    for (auto& p : model.atoms) {
        auto& atom = p.second;
        if (atom.entity_id < 0 || atom.entity_id > entities.instanced_entities.size()) continue;

        auto color = getColorFromSymbol(atom.symbol);
        if (!strcmp(atom.res_name, "HOH")) {
            if (!settings.draw_water_atoms) continue;
            color.w = settings.water_atoms_alpha;
        }
        else if (atom.is_heterogen) {
            if (!settings.draw_hetero_atoms) continue;
            color.w = settings.hetero_atoms_alpha;
        }
        else {
            if (!settings.draw_residue_atoms) continue;
            color.w = settings.residue_atoms_alpha;
        }

        entities.instanced_entities[atom.entity_id].albedo = color;
    }

    if (settings.draw_residue_ribbons) {
        for (auto& p : model.chains) {
            auto& chain = p.second;

            auto chain_color = glm::vec4(randomSaturatedColor((unsigned int)chain.chain_id), 1.0);
            for (auto& r : chain.residues) {
                if (r->entity_id < 0 || r->entity_id > entities.mesh_entities.size()) continue;

                auto& m_e = entities.mesh_entities[r->entity_id];
                switch (settings.residue_color_mode)
                {
                case PdbResidueColorMode::SECONDARY:
                {
                    switch (r->type) {
                    case PdbResidueType::COIL:
                    {
                        m_e.albedo = coil_color;
                        break;
                    }
                    case PdbResidueType::HELIX:
                    {
                        m_e.albedo = helix_color;
                        break;
                    }
                    case PdbResidueType::STRAND:
                    {
                        m_e.albedo = strand_color;
                        break;
                    }
                    }
                    break;
                }
                case PdbResidueColorMode::CHAIN:
                {
                    m_e.albedo = chain_color;
                    break;
                }
                case PdbResidueColorMode::AMINO_ACID:
                {
                    m_e.albedo = getColorFromResidueName(r->res_name);
                    break;
                }
                default:
                {
                    fprintf(stderr, "unknown residue color mode\n");
                    break;
                }
                }
                m_e.albedo.w = settings.residue_ribbons_alpha;
            }
        }
    }
}

void createEntitiesFromPdbModel(Entities &entities, PdbModel &model, PdbDrawSettings &settings, Camera &camera) {
    static const float relative_cylinder_size = 0.10;
    static const float relative_sphere_size   = 0.15;
    float avg_bond_length = 0.0;

    // We need to calculate average bond length before hand to determine size, 
    // perhaps this could be done by modifying size matrix
    for(const auto &p : model.connections) {
        auto &b = p.second;

        auto atom_1_lu = model.atoms.find(b.atom_1_id);
        if(atom_1_lu == model.atoms.end()) continue;
        auto atom_1 = atom_1_lu->second;
        auto atom_2_lu = model.atoms.find(b.atom_2_id);
        if(atom_2_lu == model.atoms.end()) continue;
        auto atom_2 = atom_2_lu->second;

        avg_bond_length += glm::length(atom_2.position - atom_1.position);
    }
    avg_bond_length /= model.connections.size();
    float calc_cylinder_r = glm::max(avg_bond_length*relative_cylinder_size, 0.1f);
    float calc_sphere_r   = glm::max(avg_bond_length*relative_sphere_size, 0.2f);

    for(auto &p : model.connections) {
        auto &b = p.second;
        auto &type = b.type;

        auto atom_1_lu = model.atoms.find(b.atom_1_id);
        if(atom_1_lu == model.atoms.end()) continue;
        auto atom_1 = atom_1_lu->second;
        auto atom_2_lu = model.atoms.find(b.atom_2_id);
        if(atom_2_lu == model.atoms.end()) continue;
        auto atom_2 = atom_2_lu->second;

        auto color_1 = getColorFromSymbol(atom_1.symbol);
        auto color_2 = getColorFromSymbol(atom_2.symbol);

        bool is_heterogen_bond = atom_1.is_heterogen && atom_2.is_heterogen;

        if (is_heterogen_bond) {
            if (!settings.draw_hetero_atoms) continue;
            color_1.w = settings.hetero_atoms_alpha;
            color_2.w = settings.hetero_atoms_alpha;
        }
        else {
            if (!settings.draw_residue_atoms) continue;
            color_1.w = settings.residue_atoms_alpha;
            color_2.w = settings.residue_atoms_alpha;
        }

        switch(type){
            case PdbConnectionType::SINGLE:
            case PdbConnectionType::SINGLE_REDUNDANT:
            {
                b.entity_id = createSingleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            case PdbConnectionType::DOUBLE:
            case PdbConnectionType::DOUBLE_REDUNDANT:
            {

                b.entity_id = createDoubleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            case PdbConnectionType::TRIPLE:
            case PdbConnectionType::TRIPLE_REDUNDANT:
            {
                b.entity_id = createTripleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            default:
                fprintf(stderr, "Unhandled bond type %u\n", type);
                break;
        }
    }
    auto center = glm::vec3(0.0);
    for(auto &p : model.atoms) {
        auto &atom = p.second;
        center += atom.position;

        auto color = getColorFromSymbol(atom.symbol);
        if (!strcmp(atom.res_name, "HOH")) {
            if (!settings.draw_water_atoms) continue;
            color.w = settings.water_atoms_alpha;
        } else if(atom.is_heterogen) {
            if (!settings.draw_hetero_atoms) continue;
            color.w = settings.hetero_atoms_alpha;
        } else {
            if (!settings.draw_residue_atoms) continue;
            color.w = settings.residue_atoms_alpha;
        }

        atom.entity_id = createAtomEntity(atom.position, color, entities, calc_sphere_r);
    }

    // There is probably a faster algorithm for determining a bounding sphere but for now
    // just loop to find center then loop to find max distance from center
    center /= model.atoms.size();
    float max_distance = 0.0;
    for(const auto &p : model.atoms) {
        auto &a = p.second;
        auto dist = glm::length(a.position - center);
        if(dist > max_distance) {
            max_distance = dist;
        }
    }

    camera.selected = true;
    camera.selection_position = center;
    camera.selection_radius = max_distance;

    // Create entities for ribbon meshes which have been precalculated
    if (settings.draw_residue_ribbons) {
        for (auto& p : model.chains) {
            auto& chain = p.second;

            auto chain_color = glm::vec4(randomSaturatedColor((unsigned int)chain.chain_id), 1.0);
            for (auto& r : chain.residues) {
                if (!r->mesh_generated) continue;

                r->entity_id = entities.mesh_entities.size();
                auto &m_e = entities.mesh_entities.emplace_back();

                m_e.mesh = &r->mesh;
                switch (settings.residue_color_mode)
                {
                    case PdbResidueColorMode::SECONDARY:
                    {
                        switch (r->type) {
                            case PdbResidueType::COIL:
                            {
                                m_e.albedo = coil_color;
                                break;
                            }
                            case PdbResidueType::HELIX:
                            {
                                m_e.albedo = helix_color;
                                break;
                            }
                            case PdbResidueType::STRAND:
                            {
                                m_e.albedo = strand_color;
                                break;
                            }
                            default:
                            {
                                fprintf(stderr, "Unknown residue type when coloring\n");
                                break;
                            }
                        }
                        break;
                    }
                    case PdbResidueColorMode::CHAIN:
                    {
                        m_e.albedo = chain_color;
                        break;
                    }
                    case PdbResidueColorMode::AMINO_ACID:
                    {
                        m_e.albedo = getColorFromResidueName(r->res_name);
                        break;
                    }
                    default:
                    {
                        fprintf(stderr, "unknown residue color mode\n");
                        break;
                    }
                }
                m_e.albedo.w = settings.residue_ribbons_alpha;
            }
        }
    }
    model.renderable = true;
}

void PdbFile::clear() {
    models.clear();
}
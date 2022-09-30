#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <set>

#include <glm/glm.hpp>
#include "controls.hpp"
#include "glm/detail/func_geometric.hpp"
#include "glm/gtc/quaternion.hpp"

#include "globals.hpp"
#include "assets.hpp"
#include "loader.hpp"
#include "entities.hpp"
#include "graphics.hpp"
#include "utilities.hpp"

// Mostly for making presentation
constexpr bool draw_water = false;
constexpr bool draw_residue_atoms = true;
constexpr float residue_atom_alpha = 0.15;
constexpr bool draw_chains = true;
constexpr float chain_alpha = 1.0;

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

static const float sphere_r = 0.7;
static const float cylinder_r = 0.2;
static const float bond_gap_r = 0.1;

glm::vec4 getColorFromSymbol(std::string symbol) {
    auto color_1_lu = symbol_to_color_lut.find(std::string(symbol));
    if(color_1_lu == symbol_to_color_lut.end()) return glm::vec4(1.0);
    else                                        return glm::vec4(color_1_lu->second / 255.f, 1.0);
}

void createAtomEntity(const glm::vec3 &pos, const glm::vec4 &col, Entities &entities, float radius=sphere_r) {
    auto &m_e = entities.instanced_entities.emplace_back(entities.instanced_entities.size());

    m_e.albedo = col;
    m_e.instance_mesh = (int)InstanceNames::SPHERE;
    m_e.scale = glm::mat3(radius);
    m_e.position = pos;
}

void createCylinderEntity(const glm::vec3 &pos_1, const glm::vec3 &pos_2, const glm::vec4 &col, Entities &entities, float radius=cylinder_r) {
    auto delta = pos_2 - pos_1;
    auto distance = glm::length(delta); 

    auto &e = entities.instanced_entities.emplace_back(entities.instanced_entities.size());
    e.instance_mesh = (int)InstanceNames::CYLINDER;
    scaleMat3(e.scale, glm::vec3(distance,  radius, radius));

    e.position = pos_1;
    e.albedo = col;
    e.rotation = quatAlignAxisToDirection(glm::vec3(1,0,0),  delta);
}

void createSingleBondEntities(const glm::vec3 &pos_1, const glm::vec4 &col_1, const glm::vec3 &pos_2, const glm::vec4 &col_2, Entities &entities, float radius=cylinder_r) {
    auto delta = pos_2 - pos_1;
    auto distance = glm::length(delta); 

    for(int i = 0; i < 2; ++i) {
        auto &e = entities.instanced_entities.emplace_back(entities.instanced_entities.size());
        e.instance_mesh = (int)InstanceNames::CYLINDER;
        scaleMat3(e.scale, glm::vec3(distance/2.0,  radius, radius));

        if(i == 0) {
            e.position = pos_1;
            e.albedo = col_1;
            e.rotation = quatAlignAxisToDirection(glm::vec3(1,0,0),  delta);
        } else {
            e.position = pos_2;
            e.albedo = col_2;
            e.rotation = quatAlignAxisToDirection(glm::vec3(1,0,0), -delta);
        }
    }
}

void createDebugCartesian(const glm::vec3 &p, const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c, Entities &entities, float r=cylinder_r) {
    createAtomEntity(p, glm::vec4(1.0,1.0,1.0,0.5), entities, r);
    createCylinderEntity(p, p+a, glm::vec4(1,0,0,0.5), entities,r);
    createCylinderEntity(p, p+b, glm::vec4(0,1,0,0.5), entities,r);
    createCylinderEntity(p, p+c, glm::vec4(0,0,1,0.5), entities,r);
}

void createDoubleBondEntities(const glm::vec3 &pos_1, const glm::vec4 &col_1, const glm::vec3 &pos_2, const glm::vec4 &col_2, Entities &entities, float radius=cylinder_r) {
    auto v_offset = anyPerpendicular(pos_2 - pos_1);
    for(int j = 0; j < 2; ++j){
        auto offset_1 = pos_1 + bond_gap_r*v_offset;
        auto offset_2 = pos_2 + bond_gap_r*v_offset;
        createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities, radius/2.0);
        v_offset *= -1.0;
    }
}

void createTripleBondEntities(const glm::vec3 &pos_1, const glm::vec4 &col_1, const glm::vec3 &pos_2, const glm::vec4 &col_2, Entities &entities, float radius=cylinder_r) {
    auto v_offset = anyPerpendicular(pos_2 - pos_1);
    for(int j = 0; j < 3; ++j){
        auto offset_1 = pos_1 + bond_gap_r*v_offset;
        auto offset_2 = pos_2 + bond_gap_r*v_offset;
        createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities, radius/2.0);
        v_offset = v_offset*glm::angleAxis((float)(2.0/3.0 * PI), glm::normalize(pos_2 - pos_1));
    }
}

void createEntitiesFromMolFile(Entities &entities, MolFile &data, Camera &camera){
#ifdef USE_ASSIMP
    loadMeshWithAssimp(entities.instanced_meshes.emplace_back(), "data/models/sphere.obj");
    loadMeshWithAssimp(entities.instanced_meshes.emplace_back(), "data/models/cylinder.obj");
#else
    readMeshFile(entities.instanced_meshes.emplace_back(), "data/models/sphere.mesh");
    readMeshFile(entities.instanced_meshes.emplace_back(), "data/models/cylinder.mesh");
#endif
    static const float relative_cylinder_size = 0.1;
    static const float relative_sphere_size   = 0.26;

    float avg_bond_length = 0.0;
    for(int i = 0; i < data.num_bonds; ++i){
        auto &b = data.bonds[i];
        auto type = b.type;
        auto &atom_1 = data.atoms[b.atom_1_index];
        auto &atom_2 = data.atoms[b.atom_2_index];
        avg_bond_length += glm::length(atom_2.position - atom_1.position);
    }
    avg_bond_length /= data.num_bonds;
    float calc_cylinder_r = avg_bond_length*relative_cylinder_size;
    float calc_sphere_r   = avg_bond_length*relative_sphere_size;

    auto center = glm::vec3(0.0);
    entities.instanced_entities.reserve(entities.instanced_entities.size() + data.num_atoms);
    for(int i = 0; i < data.num_atoms; ++i){
        auto &a = data.atoms[i];
        auto color = getColorFromSymbol(a.symbol);
        createAtomEntity(a.position, color, entities, calc_sphere_r);
        center += a.position;
    }
    center /= data.num_atoms;

    float max_distance = 0.0;
    for (int i = 0; i < data.num_atoms; ++i) {
        auto& a = data.atoms[i];
        auto dist = glm::length(a.position - center);
        if (dist > max_distance) {
            max_distance = dist;
        }
    }
    camera.target = center;
    camera.position = camera.target + glm::normalize(camera.position - camera.target) * max_distance * 2.2f;
    updateCameraView(camera);

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
                    model = &data.models.try_emplace(serial).first->second;
                }
                
                printf("MODEL %d\n", model->serial);
                model->serial = serial;
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

                // @note this transforms symbol to match symbol lut, better to modify lut to match pdb
                sscanf(&line[76], "%2s", atom.symbol);
                atom.symbol[1] = tolower(atom.symbol[1]);
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
                        printf("DICT CONECT: %d --> %d of type %d\n", atom_1_id, atom_2_id, ap.second.type);
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
void createPolypeptideEntity(Entities &entities, std::vector<PeptidePlane> &planes) {
    constexpr float r = 0.25;
    // @todo these should probably 
    constexpr int num_splines = 16, num_points_per_spline = 16;
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

    // [0,1] point to exavluate to get arrow head beginning 
    constexpr float arrow_length = 0.7;
    glm::vec2 arrow_profile[num_splines];
    glm::vec2 arrow_profile_normals[num_splines];
    createRectangleProfile(num_splines, arrow_profile, 10.f*r, 1.5f*r);
    createRectangleProfileNormals(num_splines, arrow_profile_normals, 10.f*r, 1.5f*r);
    glm::vec2 arrow_tip_profile[num_splines];
    glm::vec2 arrow_tip_profile_normals[num_splines];
    createRectangleProfile(num_splines, arrow_tip_profile, 0.f*r, 1.f*r);
    createRectangleProfileNormals(num_splines, arrow_tip_profile_normals, 0.f*r, 1.f*r);

    // @todo Handle small number of planes
    if(planes.size() < 3) return;
    
    // Points describing spline tube's surface
    glm::vec3 spline_tube [num_points_per_spline][num_splines];
    glm::vec3 normals_tube[num_points_per_spline][num_splines];

    const int num_planes = planes.size();

    glm::vec3 previous_normal = planes[0].normal;
    auto helix_color = glm::vec4(glm::vec3(255, 105, 180) / 255.f, chain_alpha);
    auto strand_color = glm::vec4(glm::vec3(255, 211, 0) / 255.f, chain_alpha);
    auto coil_color = glm::vec4(glm::vec3(220, 220, 220) / 255.f, chain_alpha);
    for(int i = 0; i < num_planes - 1; i++) {
        auto &m_e = entities.mesh_entities.emplace_back(entities.mesh_entities.size());

        // @note middle peptide or similar might reduce floating point errors
        // need to offset all vertex positions then
        //m_e.position = planes[i].position;

        auto &mesh = m_e.mesh;
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
        switch (planes[i].residue_1->type) {
            case PdbResidueType::COIL:
                {
                    m_e.albedo = coil_color;
                    pf1  = circle_profile;
                    pfn1 = circle_profile_normals;
                    break;
                }
            case PdbResidueType::HELIX:
                {
                    m_e.albedo = helix_color;
                    pf1  = ribbon_profile;
                    pfn1 = ribbon_profile_normals;
                    break;
                }
            case PdbResidueType::STRAND:
                {
                    m_e.albedo = strand_color;
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
            auto arrow_base_plane = createPartialHermiteSplineNormalsBetweenProfiles(num_points_per_spline, num_splines, 
                    rectangle_profile, rectangle_profile_normals, rectangle_profile, rectangle_profile_normals, 
                    planes[p1_i], planes[i], planes[i+1], planes[p4_i], 
                    &spline_tube[0][0], 
                    &normals_tube[0][0], previous_normal, arrow_length);
            createClosedSurfaceFromSplinesNormals(num_splines, num_points_per_spline, &spline_tube[0][0], &normals_tube[0][0], mesh);

            glm::vec3 projected_profile[num_splines];
            projectPointsOnPlane(num_splines, arrow_base_plane.position, arrow_base_plane.normal, arrow_base_plane.right, 
                    arrow_profile, projected_profile);
            createClosedFacedFromProfile(arrow_base_plane.position, num_splines, projected_profile, mesh, true);

            if(i == num_planes - 2) {
                // hack to make longer arrows at end points
                planes[i+1].position += arrow_base_plane.forward;
                createHermiteSplineNormalsBetweenProfiles(2, num_splines, 
                        arrow_profile, arrow_profile_normals, arrow_tip_profile, arrow_tip_profile_normals, 
                        arrow_base_plane, arrow_base_plane, planes[i+1], planes[i+1], 
                        &spline_tube[0][0], 
                        &normals_tube[0][0], previous_normal);
            } else {
                createHermiteSplineNormalsBetweenProfiles(2, num_splines, 
                        arrow_profile, arrow_profile_normals, circle_profile, circle_profile_normals, 
                        arrow_base_plane, arrow_base_plane, planes[i+1], planes[i+1], 
                        &spline_tube[0][0], 
                        &normals_tube[0][0], previous_normal);
            }
            createClosedSurfaceFromSplinesNormals(num_splines, 2, &spline_tube[0][0], &normals_tube[0][0], mesh);

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

            createHermiteSplineNormalsBetweenProfiles(num_points_per_spline, num_splines, pf1, pfn1, pf2, pfn2, 
                    planes[p1_i], planes[i], planes[i+1], planes[p4_i], &spline_tube[0][0], 
                    &normals_tube[0][0], previous_normal);
            createClosedSurfaceFromSplinesNormals(num_splines, num_points_per_spline, &spline_tube[0][0], &normals_tube[0][0], mesh);

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
    }
}

void createEntitiesFromPdbFile(Entities &entities, PdbFile &data, Camera &camera){
#ifdef USE_ASSIMP
    loadMeshWithAssimp(entities.instanced_meshes.emplace_back(), "data/models/sphere.obj");
    loadMeshWithAssimp(entities.instanced_meshes.emplace_back(), "data/models/cylinder.obj");
#else
    readMeshFile(entities.instanced_meshes.emplace_back(), "data/models/sphere.mesh");
    readMeshFile(entities.instanced_meshes.emplace_back(), "data/models/cylinder.mesh");
#endif

    // @note Assume first model is the correct one
    if (data.models.size() == 0) return;
    auto &model = data.models[0];

    static const float relative_cylinder_size = 0.10;
    static const float relative_sphere_size   = 0.15;
    float avg_bond_length = 0.0;
    for(const auto &p : model.connections) {
        auto &b = p.second;

        auto atom_1_lu = model.atoms.find(b.atom_1_id);
        if(atom_1_lu == model.atoms.end()) continue;
        auto atom_1 = atom_1_lu->second;
        auto atom_2_lu = model.atoms.find(b.atom_2_id);
        if(atom_2_lu == model.atoms.end()) continue;
        auto atom_2 = atom_2_lu->second;

        //if(!draw_residue_atoms && (!atom_1.is_heterogen || !atom_2.is_heterogen)) continue;

        avg_bond_length += glm::length(atom_2.position - atom_1.position);
    }
    avg_bond_length /= model.connections.size();
    float calc_cylinder_r = glm::max(avg_bond_length*relative_cylinder_size, 0.1f);
    float calc_sphere_r   = glm::max(avg_bond_length*relative_sphere_size, 0.2f);

    printf("Calc sphere %f, Calc cylinder %f\n", calc_sphere_r, calc_cylinder_r);

    // Draw hetero atoms and the bonds between them (including non hetero edge atoms)
    // @speed
    std::set<int> encountered_hetatms;
    for(const auto &p : model.connections) {
        auto &b = p.second;
        auto type = b.type;

        auto atom_1_lu = model.atoms.find(b.atom_1_id);
        if(atom_1_lu == model.atoms.end()) continue;
        auto atom_1 = atom_1_lu->second;
        auto atom_2_lu = model.atoms.find(b.atom_2_id);
        if(atom_2_lu == model.atoms.end()) continue;
        auto atom_2 = atom_2_lu->second;

        //printf("Hetero Bond %d --> %d with type %u\n", b.atom_1_id, b.atom_2_id, type);
        auto color_1 = getColorFromSymbol(atom_1.symbol);
        auto color_2 = getColorFromSymbol(atom_2.symbol);

        // Add transparency to connections between non-hetero atoms
        if(!atom_1.is_heterogen || !atom_2.is_heterogen) {
            if(!draw_residue_atoms) continue;
            color_1.w = residue_atom_alpha;
            color_2.w = residue_atom_alpha;
        } else {
            encountered_hetatms.insert(b.atom_1_id);
            encountered_hetatms.insert(b.atom_2_id);
        }

        switch(type){
            case PdbConnectionType::SINGLE:
            case PdbConnectionType::SINGLE_REDUNDANT:
            {
                createSingleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            case PdbConnectionType::DOUBLE:
            case PdbConnectionType::DOUBLE_REDUNDANT:
            {

                createDoubleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            case PdbConnectionType::TRIPLE:
            case PdbConnectionType::TRIPLE_REDUNDANT:
            {
                createTripleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            default:
                fprintf(stderr, "Unhandled bond type %u\n", type);
                break;
        }
    }
    auto center = glm::vec3(0.0);
    for(const auto &p : model.atoms) {
        auto &atom = p.second;
        center += atom.position;

        // If atom is not hetero draw it with transparency
        if(draw_residue_atoms && encountered_hetatms.find(atom.serial) == encountered_hetatms.end()) {
            // Water is a special case
            if(!draw_water && !strcmp(atom.res_name, "HOH")) continue;

            auto color = getColorFromSymbol(atom.symbol);
            color.w = residue_atom_alpha;
            //if(!strcmp(atom.name, "CA  ")) color.w = 1.0;
            createAtomEntity(atom.position, color, entities, calc_sphere_r);
        }

        //if(!a.is_heterogen) continue;
        //encountered_hetatms.insert(a.serial);
    }
    center /= model.atoms.size();
    float max_distance = 0.0;
    for(const auto &p : model.atoms) {
        auto &a = p.second;
        auto dist = glm::length(a.position - center);
        if(dist > max_distance) {
            max_distance = dist;
        }
    }
    camera.target = center;
    camera.position = camera.target + glm::normalize(camera.position - camera.target)*max_distance*2.2f;
    updateCameraView(camera);

    for(const auto &atom_id : encountered_hetatms) {
        auto &atom = model.atoms[atom_id];
        auto color = getColorFromSymbol(atom.symbol);
        createAtomEntity(atom.position, color, entities, calc_sphere_r);
    }

    if(!draw_chains) return;

    for(const auto &p : model.chains) {
        auto &chain = p.second;
        if(chain.residues.size() < 2) continue;

        std::vector<PeptidePlane> peptide_planes;
    
        // @todo adjust color to secondary structures
        //printf("Chain %c\n", chain.chain_id);
        int plane_i = 0;
        for(int i = 0; i < chain.residues.size() - 1; i++) {
            // Determine plane of connection between residues through peptide bond
            auto &r1 = chain.residues[i];
            auto &r2 = chain.residues[i+1];

            // @note Peptide bond forms between carboxylic group and amino group.
            // PDB should guarantee that residues' label the alpha carbon as "CA"
            // we just use the "O" of the amine to determine the perpendicular plane 
            // between adjacent amino acids of protein.
            auto lu_CA1 = r1->atom_name_id.find("CA  ");
            if(lu_CA1  == r1->atom_name_id.end()) continue;
            auto lu_CA2 = r2->atom_name_id.find("CA  ");
            if(lu_CA2  == r2->atom_name_id.end()) continue;
            //auto lu_O1  = r1->atom_name_id.find(" O  ");
            //if(lu_O1   == r1->atom_name_id.end()) continue;
            
            auto CA1 = model.atoms[lu_CA1->second].position;
            auto CA2 = model.atoms[lu_CA2->second].position;
            //auto O1  = model.atoms[lu_O1 ->second].position;

            auto &plane = peptide_planes.emplace_back();
            plane.residue_1 = r1;
            plane.residue_2 = r2;

            plane.position = (CA1 + CA2) / 2.f;

            plane.forward = glm::normalize(CA2 - CA1);
            //plane.normal  = glm::normalize(glm::cross(plane.forward, O1 - CA1));
            //// To ensure perpendicularity calculate right
            //plane.right   = glm::normalize(glm::cross(plane.normal, plane.forward));

            plane.normal = glm::vec3(0);

            //printf("Peptide Plane %d --> %d\n", i, i+1);
        }
        /*if(peptide_planes.size() > 0) {
            auto plane = peptide_planes[peptide_planes.size() - 1];
            auto lu_CA2 = plane.residue_2->atom_name_id.find(" CA ");
            plane.position = model.atoms[lu_CA2->second].position;
            peptide_planes.emplace_back(plane);
        }*/
        // https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Representations
        // cardinal spline calculation
        constexpr float c = 0.25;
        for(int i = 1; i < peptide_planes.size() - 1; i++) {
            auto &p = peptide_planes;
            p[i].forward = (1-c) * (p[i+1].position - p[i-1].position);
        }
        for(int i = 0; i < peptide_planes.size() - 1; i++) {
            auto &p = peptide_planes;
            auto n0 = -6.f*p[i].position - 4.f*p[i].forward + 6.f*p[i+1].position - 2.f*p[i+1].forward;
            n0 = perpendicularComponent(n0, p[i].forward);
            p[i].normal += n0 / 2.f;
            auto n1 =  6.f*p[i].position + 2.f*p[i].forward - 6.f*p[i+1].position + 4.f*p[i+1].forward;
            n1 = perpendicularComponent(n1, p[i+1].forward);
            p[i+1].normal += n1 / 2.f;
        }
        for(int i = 0; i < peptide_planes.size(); i++) {
            auto &p = peptide_planes;
            p[i].normal = glm::normalize(p[i].normal);
            if(i > 0){
                if(glm::dot(p[i].normal, p[i-1].normal) < 0) {
                    p[i].normal *= -1;
                }
            }
            p[i].right = glm::normalize(glm::cross(p[i].forward, p[i].normal));
        }
        if(peptide_planes.size() > 1) {
            auto &p = peptide_planes;
            p[0].normal = glm::normalize(glm::cross(p[0].forward, p[0].right));
            p[peptide_planes.size()-1].normal = glm::normalize(glm::cross(p[peptide_planes.size()-1].forward, p[peptide_planes.size()-1].right));
        }

        // Flip peptides that are >180 from previous
        for(int i = 1; i < peptide_planes.size(); i++) {
            if(glm::dot(peptide_planes[i].normal, peptide_planes[i - 1].normal) < 0) {
                peptide_planes[i].normal *= -1;
                peptide_planes[i].right  *= -1;
                peptide_planes[i].flipped = true;
            }
        }

        //for(int i = 1; i < peptide_planes.size(); i++) {
            //createDebugCartesian(plane.position, plane.normal)
            //createAtomEntity(peptide_planes[i].position, glm::vec4(1.0, 0.0, 0.0, 1.0), entities, calc_sphere_r);
            //auto fwd = glm::normalize(peptide_planes[i].forward)*1.5f;
            //createSingleBondEntities(peptide_planes[i].position - fwd*0.5f, glm::vec4(1,0,0,1), peptide_planes[i].position+fwd*0.5f, glm::vec4(1,0,0,1), entities,calc_cylinder_r/2.0f);
        //}

        createPolypeptideEntity(entities, peptide_planes);
    }
}

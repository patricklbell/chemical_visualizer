#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <set>

#include <glm/glm.hpp>
#include "glm/detail/func_geometric.hpp"
#include "glm/gtc/quaternion.hpp"

#include "assets.hpp"
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

static const float sphere_r = 0.7;
static const float cylinder_r = 0.2;
static const float bond_gap_r = 0.1;

glm::vec3 getColorFromSymbol(std::string symbol) {
    auto color_1_lu = symbol_to_color_lut.find(std::string(symbol));
    if(color_1_lu == symbol_to_color_lut.end()) return glm::normalize(glm::vec3(1.0));
    else                                        return glm::normalize(color_1_lu->second);
}



void createSingleBondEntities(const glm::vec3 &pos_1, const glm::vec3 &col_1, const glm::vec3 &pos_2, const glm::vec3 &col_2, Entities &entities, float radius=cylinder_r) {
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
    createSingleBondEntities(p, glm::vec3(1,0,0), p+a, glm::vec3(1,0,0), entities,r);
    createSingleBondEntities(p, glm::vec3(0,1,0), p+b, glm::vec3(0,1,0), entities,r);
    createSingleBondEntities(p, glm::vec3(0,0,1), p+c, glm::vec3(0,0,1), entities,r);
}

void createDoubleBondEntities(const glm::vec3 &pos_1, const glm::vec3 &col_1, const glm::vec3 &pos_2, const glm::vec3 &col_2, Entities &entities, float radius=cylinder_r) {
    auto v_offset = anyPerpendicular(pos_2 - pos_1);
    for(int j = 0; j < 2; ++j){
        auto offset_1 = pos_1 + bond_gap_r*v_offset;
        auto offset_2 = pos_2 + bond_gap_r*v_offset;
        createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities, radius);
        v_offset *= -1.0;
    }
}

void createTripleBondEntities(const glm::vec3 &pos_1, const glm::vec3 &col_1, const glm::vec3 &pos_2, const glm::vec3 &col_2, Entities &entities, float radius=cylinder_r) {
    auto v_offset = anyPerpendicular(pos_2 - pos_1);
    for(int j = 0; j < 3; ++j){
        auto offset_1 = pos_1 + bond_gap_r*v_offset;
        auto offset_2 = pos_2 + bond_gap_r*v_offset;
        createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities, radius);
        v_offset = v_offset*glm::angleAxis((float)(2.0/3.0 * PI), glm::normalize(pos_2 - pos_1));
    }
}

void createAtomEntity(const glm::vec3 &pos, const glm::vec3 &col, Entities &entities, float radius=sphere_r) {
    auto &m_e = entities.instanced_entities.emplace_back(entities.instanced_entities.size());

    m_e.albedo = col;
    m_e.instance_mesh = (int)InstanceNames::SPHERE;
    m_e.scale = glm::mat3(radius);
    m_e.position = pos;
}

glm::vec3 createEntitiesFromMolFile(Entities &entities, MolFile &data){
    // @todo make these mesh files to load quicker
    loadMeshWithAssimp(entities.instanced_meshes.emplace_back(), "data/models/sphere.obj");
    loadMeshWithAssimp(entities.instanced_meshes.emplace_back(), "data/models/cylinder.obj");

    static const float relative_cylinder_size = 0.05;
    static const float relative_sphere_size   = 0.3;

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
    return center;
}

long hashBondPair(int id_1, int id_2) {
    if(id_2 < id_1) return (long)id_1 ^ ((long)id_2 << 16);
    else            return (long)id_2 ^ ((long)id_1 << 16);
}

void createResidueId(int seq_num, char res_name[4], char chain_id, char i_code, char result[10]) {
    sprintf(result, "%3s%c%4d%c", res_name, chain_id, seq_num, i_code); 
}

void loadPdbFile(PdbFile &data, std::string path){
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

    auto &model = data.models.emplace_back();
    model.serial = 0;
    // Handle model command being missing
    bool first_model_flag = true;

    while(1) {
        if(fgets(line, 1024, f) == NULL) break;

        if(strlen(line) >= 6) {
            substrString(line, 0, 5, record_name);

            if(!strcmp(record_name, "MODEL ")) {
                if(first_model_flag) {
                    first_model_flag = false;
                } else {
                    model = data.models.emplace_back();
                }
                sscanf(line, "MODEL %d", &model.serial);
                printf("MODEL %d\n", model.serial);
            } else if(!strcmp(record_name, "HETATM") || !strcmp(record_name, "ATOM  ")) {
                int serial;
                substrSscanf(line, 6, 10, " %d", &serial);

                // @note this overwrites duplicate ids
                PdbAtom &atom = model.atoms[serial];
                atom.serial = serial;
                atom.is_heterogen = !strcmp(record_name, "HETATM");

                substrString(line, 12, 15,  atom.name);

                substrString(line, 17, 19,  atom.res_name);
                substrChar  (line, 21,     &atom.chain_id);
                substrInt   (line, 22, 25, &atom.res_seq);

                auto lu = model.residues.find(atom.res_seq);
                // If residue doesn't exist create it
                if(lu == model.residues.end()) {
                    auto &residue = model.residues.try_emplace(atom.res_seq).first->second;

                    memcpy(residue.res_name, atom.res_name, 4);
                    residue.chain_id = atom.chain_id;
                    residue.i_code   = atom.i_code;
                    residue.res_seq  = atom.res_seq;

                    // Add new residue to chain
                    auto lu = model.chains.find(atom.chain_id);
                    // If chain doesn't exist create it
                    // @note doesn't account for empty chain_id i.e ' '
                    if(lu == model.chains.end()) {
                        auto &chain = model.chains.try_emplace(atom.chain_id).first->second;
                        chain.chain_id = atom.chain_id;
                        // @note pointers to values in unordered_map are stable 
                        chain.residues.push_back(&residue);

                        printf("CHAIN : %c\n", chain.chain_id); // @debug
                    } else {
                        auto &chain = lu->second;
                        chain.residues.push_back(&residue);
                    }

                    residue.atom_ids.push_back(atom.serial);
                    residue.atom_name_id[atom.name] = atom.serial;

                    printf("RESDUE: res_name %s chain_id %c i_code %c res_seq %d\n", 
                            residue.res_name, residue.chain_id, residue.i_code, residue.res_seq); // @debug
                } else {
                    auto &residue = lu->second;
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
                    auto lu = model.connections.find(hash);
                    if(lu != model.connections.end()) {
                        auto &b = lu->second;
                        b.type = PdbConnectionType((unsigned int)b.type + 1);
                    } else {
                        auto b = &model.connections[hash];
                        b->atom_1_id = atom_id;
                        b->atom_2_id = connect_ids[i];
                        b->type = PdbConnectionType::SINGLE;
                    }
                }
                printf("\n");
            } else if(!strcmp(record_name, "HELIX ")) { 
                int serial;
                substrInt(line, 7, 9, &serial);
                auto &helix = model.helices.try_emplace(serial).first->second;
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
                auto lu = model.sheets.find(sheet_id);
                PdbSheet *sheet;
                if(lu == model.sheets.end()) {
                    // Create and modify sheet only first time
                    sheet = &model.sheets.try_emplace(sheet_id).first->second;
                    memcpy(sheet->sheet_id, sheet_id, 4);
                    substrInt(   line, 14, 15, &sheet->num_strands);
                    bool is_first = true;
                } else {
                    sheet = &model.sheets[sheet_id];
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
            } else if(!strcmp(record_name, "ENDMDL")) {
                printf("ENDMDL\n");
            }
        }
    }

    // Processing after loading from file
    for(auto &model : data.models) {
        // Label residue type based on secondary structures
        for(auto &j : model.helices) {
            auto &helix = j.second;
            for(int seq_num = helix.init_seq_num; seq_num <= helix.end_seq_num; seq_num++) {
                auto lu = model.residues.find(seq_num);
                if(lu == model.residues.end()) continue;

                auto &residue = lu->second;
                residue.type = PdbResidueType::HELIX;
            }
        }
        for(auto &j : model.sheets) {
            auto &sheet = j.second;
            for(auto &k : sheet.strands) {
                auto &strand = k.second;
                for(int seq_num = strand.init_seq_num; seq_num <= strand.end_seq_num; seq_num++) {
                    auto lu = model.residues.find(seq_num);
                    if(lu == model.residues.end()) continue;

                    auto &residue = lu->second;
                    residue.type = PdbResidueType::STRAND;
                }
            }
        }
    }
}

// Mesh which blends between different secondary structures in one polypeptide chain 
void createPolypeptideEntity(Entities &entities, std::vector<PeptidePlane> &planes) {
    auto color = randomColor();

    auto &m_e = entities.mesh_entities.emplace_back(entities.mesh_entities.size());

    // @note middle peptide or similar might reduce floating point errors
    // need to offset all vertex positions then
    //m_e.position = planes[0].position;
    m_e.albedo = randomColor();

    auto &mesh = m_e.mesh;
	mesh.draw_mode = GL_TRIANGLES;
	mesh.draw_type = GL_UNSIGNED_SHORT;

    // @todo make parts of mesh different colors
    mesh.num_materials = 1;
    mesh.draw_start = (GLint*)malloc(sizeof(GLint) * mesh.num_materials);
    mesh.draw_count = (GLint*)malloc(sizeof(GLint) * mesh.num_materials);
    mesh.draw_start[0] = 0;
    mesh.draw_count[0] = 0;
    printf("Created Peptide mesh\n");

    constexpr float r = 0.3;
    // @todo these should probably 
    constexpr int num_splines = 20, num_points_per_spline = 5;
    glm::vec2 circle_profile[num_splines];
    createCircularProfile(num_splines, circle_profile, r);

    printf("Created circular profile\n");
    for(int i = 0; i < num_splines; i++) {
        auto &pf = circle_profile[i];
        printf("Vertex %3d x: %9.6f y: %9.6f\n", i, pf.x, pf.y);
    }




    // @debug -->
    glm::vec3 tube[num_splines][num_points_per_spline];
    
    glm::vec3 pf1[num_splines];
    glm::vec3 pf2[num_splines];
    glm::vec3 pf3[num_splines];
    glm::vec3 pf4[num_splines];
    //auto n = glm::normalize(glm::vec3(rand() - RAND_MAX / 2, rand() - RAND_MAX / 2, rand() - RAND_MAX / 2));
    //auto u = glm::normalize(anyPerpendicular(n));
    //auto v = glm::normalize(glm::cross(n, u));
    projectPointsOnPlane(num_splines, glm::vec3(-0.5,0,0), glm::vec3(0,0,1), glm::vec3(0,1,0), circle_profile, pf1);
    projectPointsOnPlane(num_splines, glm::vec3( 0.0,0,0), glm::vec3(0,0,1), glm::vec3(0,1,0), circle_profile, pf2);
    projectPointsOnPlane(num_splines, glm::vec3( 0.5,0,0), glm::vec3(0,0,1), glm::vec3(0,1,0), circle_profile, pf3);
    projectPointsOnPlane(num_splines, glm::vec3( 1.0,0,0), glm::vec3(0,0,1), glm::vec3(0,1,0), circle_profile, pf4);
    //projectPointsOnPlane(num_splines, glm::vec3(0,0,0), u, v, circle_profile, pf4);
    createClosedFacedFromProfile(glm::vec3(-0.5,0,0), num_splines, pf1, mesh, false);
    createClosedFacedFromProfile(glm::vec3( 0.0,0,0), num_splines, pf2, mesh, false);
    createClosedFacedFromProfile(glm::vec3( 0.5,0,0), num_splines, pf3, mesh, false);
    createClosedFacedFromProfile(glm::vec3( 1.0,0,0), num_splines, pf4, mesh, false);
    //createClosedFacedFromProfile(glm::vec3(0,0,0), num_splines, pf4, mesh);
    createDebugCartesian(glm::vec3(-0.5,0,0), 0.2f*glm::vec3(1,0,0), 0.2f*glm::vec3(0,0,1), 0.2f*glm::vec3(0,1,0), entities, 0.03);
    createDebugCartesian(glm::vec3( 0.0,0,0), 0.2f*glm::vec3(1,0,0), 0.2f*glm::vec3(0,0,1), 0.2f*glm::vec3(0,1,0), entities, 0.03);
    createDebugCartesian(glm::vec3( 0.5,0,0), 0.2f*glm::vec3(1,0,0), 0.2f*glm::vec3(0,0,1), 0.2f*glm::vec3(0,1,0), entities, 0.03);
    createDebugCartesian(glm::vec3( 1.0,0,0), 0.2f*glm::vec3(1,0,0), 0.2f*glm::vec3(0,0,1), 0.2f*glm::vec3(0,1,0), entities, 0.03);
    //createDebugCartesian(glm::vec3(0,0,0), 0.2f*n, 0.2f*u, 0.2f*v, entities, 0.03);

    for(int j = 0; j < num_splines; j++) {
        // Create spline between corresponding points of adjacent profiles
        createCubicSpline(pf1[j], pf2[j], pf3[j], pf4[j], num_points_per_spline, tube[j]);
    }

    createClosedSurfaceFromSplines(num_splines, num_points_per_spline, (glm::vec3*)tube, mesh);

    createMeshVao(mesh);
    return;
    // <-- @debug




    
    // Moving window of profiles
    glm::vec3 profiles[num_splines][4];

    // Points describing spline tube's surface
    glm::vec3 spline_tube[num_splines][num_points_per_spline];
    for(int i = 0; i < planes.size() - 2; i++) {
        // Keeps track of the moving window, this is element "0"
        int offset = i % 4;
        auto &pf1 = profiles[(offset  ) % 4];
        auto &pf2 = profiles[(offset+1) % 4];
        auto &pf3 = profiles[(offset+2) % 4];
        auto &pf4 = profiles[(offset+3) % 4];

        if(i == 0) {
            // Fake p1 to handle start point 
            auto &p2 = planes[0];
            auto &p3 = planes[1];
            auto &p4 = planes[2];
            projectPointsOnPlane(num_splines, p2.position, p2.right, p2.forward, circle_profile, profiles[1]);
            projectPointsOnPlane(num_splines, p3.position, p3.right, p3.forward, circle_profile, profiles[2]);
            projectPointsOnPlane(num_splines, p4.position, p4.right, p4.forward, circle_profile, profiles[3]);

            // Kinda hacky way of making natural spline beginning by preserving slope 
            for(int j = 0; j < num_splines; j++) {
                pf1[j] = -(pf3[j] - pf2[j]) + pf2[j];
            }
        } else if (i == planes.size() - 2){
            // Same hack for endpoint to preserve slope
            for(int j = 0; j < num_splines; j++) {
                pf4[j] = (pf3[j] - pf2[j]) + pf3[j];
            }
        } else {
            auto &p4 = planes[i + 3];
            // Update the "4th" profile with 4th plane/final control point
            projectPointsOnPlane(num_splines, p4.position, p4.right, p4.forward, circle_profile, pf4);
        }
        glm::vec3 test_tube[num_splines][2];
        for(int j = 0; j < num_splines; j++) {
            test_tube[j][0] = pf2[j];
            test_tube[j][1] = pf3[j];
        }
        createClosedSurfaceFromSplines(num_splines, 2, (glm::vec3*)test_tube, mesh);

        //for(int j = 0; j < num_splines; j++) {
        //    // Create spline between corresponding points of adjacent profiles
        //    createCubicSpline(pf1[j], pf2[j], pf3[j], pf4[j], num_points_per_spline, spline_tube[j]);
        //}

        //createClosedSurfaceFromSplines(num_splines, num_points_per_spline, (glm::vec3*)spline_tube, mesh);
    }
    createMeshVao(mesh);
}

glm::vec3 createEntitiesFromPdbFile(Entities &entities, PdbFile &data){
    // @todo make these mesh files to load quicker
    loadMeshWithAssimp(entities.instanced_meshes.emplace_back(), "data/models/sphere.obj");
    loadMeshWithAssimp(entities.instanced_meshes.emplace_back(), "data/models/cylinder.obj");

    auto center = glm::vec3(0.0);

    // @note Assume first model is the correct one
    if (data.models.size() == 0) return glm::vec3(0,0,0);
    auto &model = data.models[0];

    static const float relative_cylinder_size = 0.05;
    static const float relative_sphere_size   = 0.3;
    float avg_bond_length = 0.0;
    for(const auto &p : model.connections) {
        auto &b = p.second;

        auto atom_1_lu = model.atoms.find(b.atom_1_id);
        if(atom_1_lu == model.atoms.end()) continue;
        auto atom_1 = atom_1_lu->second;
        auto atom_2_lu = model.atoms.find(b.atom_2_id);
        if(atom_2_lu == model.atoms.end()) continue;
        auto atom_2 = atom_2_lu->second;

        if(!atom_1.is_heterogen && !atom_2.is_heterogen) continue;

        avg_bond_length += glm::length(atom_2.position - atom_1.position);
    }
    avg_bond_length /= model.connections.size();
    float calc_cylinder_r = avg_bond_length*relative_cylinder_size;
    float calc_sphere_r   = avg_bond_length*relative_sphere_size;

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

        if(!atom_1.is_heterogen && !atom_2.is_heterogen) continue;
        encountered_hetatms.insert(b.atom_1_id);
        encountered_hetatms.insert(b.atom_2_id);

        printf("Hetero Bond %d --> %d with type %u\n", b.atom_1_id, b.atom_2_id, type);
        auto color_1 = getColorFromSymbol(atom_1.symbol);
        auto color_2 = getColorFromSymbol(atom_2.symbol);

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
    for(const auto &p : model.atoms) {
        auto &a = p.second;
        center += a.position;

        if(!a.is_heterogen) continue;
        encountered_hetatms.insert(a.serial);
    }
    center /= model.atoms.size();

    for(const auto &atom_id : encountered_hetatms) {
        auto &atom = model.atoms[atom_id];
        auto color = getColorFromSymbol(atom.symbol);
        createAtomEntity(atom.position, color, entities, calc_sphere_r);
    }

    // @todo create list of planes of helix from atoms if they are not 1/2 not "CA" or 1 not "O"
    for(const auto &p : model.chains) {
        auto chain = p.second;
        if(chain.residues.size() < 2) continue;

        std::vector<PeptidePlane> peptide_planes;
    
        // @todo adjust color to secondary structures
        printf("Chain %d\n", chain.chain_id);
        for(int i = 0; i < chain.residues.size() - 1; i++) {
            auto &plane = peptide_planes.emplace_back();

            // Determine plane of connection between residues through peptide bond
            plane.residue_1 = chain.residues[i];
            auto &r1 = plane.residue_1;
            plane.residue_2 = chain.residues[i+1];
            auto &r2 = plane.residue_2;

            // @note Peptide bond forms between carboxylic group and amino group.
            // PDB should guarantee that residues' O, C in peptide bond are named
            // "O", "CA" so we can determine plane of residue from relationship
            // between peptides of adjacent amino acids.
            auto lu_CA1 = r1->atom_name_id.find(" CA ");
            if(lu_CA1  == r1->atom_name_id.end()) continue;
            auto lu_CA2 = r2->atom_name_id.find(" CA ");
            if(lu_CA2  == r2->atom_name_id.end()) continue;
            auto lu_O1  = r1->atom_name_id.find(" O  ");
            if(lu_O1   == r1->atom_name_id.end()) continue;
            
            auto CA1 = model.atoms[lu_CA1->second].position;
            auto CA2 = model.atoms[lu_CA2->second].position;
            auto O1  = model.atoms[lu_O1 ->second].position;

            // Average positions of peptide bonds is approximate center of residue
            plane.position = (CA1 + CA2) / 2.0f;

            plane.forward = glm::normalize(CA2 - CA1);
            plane.normal  = glm::normalize(glm::cross(plane.forward, O1 - CA1));
            // To ensure perpendicularity calculate right
            plane.right   = glm::normalize(glm::cross(plane.normal, plane.forward));
            printf("Peptide Plane %d --> %d\n", i, i+1);
        }
        // Flip peptides that are >180 from previous
        for(int i = 1; i < peptide_planes.size(); i++) {
            if(glm::dot(peptide_planes[i].normal, peptide_planes[i - 1].normal) < 0) {
                peptide_planes[i].normal *= -1;
                peptide_planes[i].right  *= -1;
                peptide_planes[i].flipped = true;
            }
        }

        createPolypeptideEntity(entities, peptide_planes);
        return glm::vec3(0,0,0);
    }

    return center;
}

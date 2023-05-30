#include <loader.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <istream>
#include <unordered_map>
#include <set>

#include <glm/glm.hpp>
#include <glm/detail/func_geometric.hpp>
#include <glm/gtc/quaternion.hpp>

#include <utilities/math.hpp>
#include <utilities/colors.hpp>
#include <utilities/strings.hpp>
#include <utilities/mesh_generation.hpp>
#include <camera/core.hpp>
#include <mesh.hpp>
#include <entities.hpp>

//
// General helper functions used by all loaders
//

static const float sphere_r   = 0.7;
static const float cylinder_r = 0.2;
static const float bond_gap_r = 0.1;

void createDebugCartesian(const glm::vec3 &p, const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c, Entities &entities, float r);
void createDebugCartesian(const Frame &f, Entities &entities, float r);

glm::vec4 getColorFromSymbol(std::string symbol) {
    auto color_1_lu = symbol_to_color_lut.find(std::string(symbol));
    if (color_1_lu == symbol_to_color_lut.end()) {
        std::cerr << "Unknown element symbol " << symbol << ".\n";
        return glm::vec4(1.0);
    } else {
        return glm::vec4(color_1_lu->second / 255.f, 1.0);
    }
}

int createAtomEntity(const glm::vec3 &pos, const glm::vec4 &col, Entities &entities, float radius = sphere_r) {
    int id    = entities.instanced_entities.size();
    auto &m_e = entities.instanced_entities.emplace_back();

    m_e.albedo   = col;
    m_e.type     = InstancedEntity::Type::SPHERE;
    m_e.scale    = glm::mat3(radius);
    m_e.position = pos;

    return id;
}

int createCylinderEntity(const glm::vec3 &pos_1, const glm::vec3 &pos_2, const glm::vec4 &col, Entities &entities,
                         float radius = cylinder_r) {
    auto delta    = pos_2 - pos_1;
    auto distance = glm::length(delta);

    int id  = entities.instanced_entities.size();
    auto &e = entities.instanced_entities.emplace_back();
    e.type  = InstancedEntity::Type::CYLINDER;
    scale_mat3(e.scale, glm::vec3(distance, radius, radius));

    e.position = pos_1;
    e.albedo   = col;
    e.rotation = rotate_to(glm::vec3(1, 0, 0), delta);

    return id;
}

int createSingleBondEntities(const glm::vec3 &pos_1, const glm::vec4 &col_1, const glm::vec3 &pos_2,
                             const glm::vec4 &col_2, Entities &entities, float radius = cylinder_r) {
    auto delta    = pos_2 - pos_1;
    auto distance = glm::length(delta);

    int id;
    for (int i = 0; i < 2; ++i) {
        id      = entities.instanced_entities.size();
        auto &e = entities.instanced_entities.emplace_back();
        e.type  = InstancedEntity::Type::CYLINDER;
        scale_mat3(e.scale, glm::vec3(distance / 2.0, radius, radius));

        if (i == 0) {
            e.position = pos_1;
            e.albedo   = col_1;
            e.rotation = rotate_to(glm::vec3(1, 0, 0), delta);
        } else {
            e.position = pos_2;
            e.albedo   = col_2;
            e.rotation = rotate_to(glm::vec3(1, 0, 0), -delta);
        }
    }
    return id;
}

void createDebugCartesian(const glm::vec3 &p, const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c,
                          Entities &entities, float r = cylinder_r) {
    createAtomEntity(p, glm::vec4(1.0, 1.0, 1.0, 0.5), entities, r);
    createCylinderEntity(p, p + a, glm::vec4(1, 0, 0, 0.5), entities, r);
    createCylinderEntity(p, p + b, glm::vec4(0, 1, 0, 0.5), entities, r);
    createCylinderEntity(p, p + c, glm::vec4(0, 0, 1, 0.5), entities, r);
}

void createDebugCartesian(const Frame &f, Entities &entities, float r = cylinder_r) {
    createDebugCartesian(f.position, f.tangent, f.normal, f.binormal, entities, r);
}

int createDoubleBondEntities(const glm::vec3 &pos_1, const glm::vec4 &col_1, const glm::vec3 &pos_2,
                             const glm::vec4 &col_2, Entities &entities, float radius = cylinder_r) {
    auto v_offset = any_perpendicular(pos_2 - pos_1);
    int id;
    for (int j = 0; j < 2; ++j) {
        auto offset_1 = pos_1 + bond_gap_r * v_offset;
        auto offset_2 = pos_2 + bond_gap_r * v_offset;
        id            = createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities, radius / 2.0);
        v_offset *= -1.0;
    }
    return id;
}

int createTripleBondEntities(const glm::vec3 &pos_1, const glm::vec4 &col_1, const glm::vec3 &pos_2,
                             const glm::vec4 &col_2, Entities &entities, float radius = cylinder_r) {
    auto v_offset = any_perpendicular(pos_2 - pos_1);
    int id;
    for (int j = 0; j < 3; ++j) {
        auto offset_1 = pos_1 + bond_gap_r * v_offset;
        auto offset_2 = pos_2 + bond_gap_r * v_offset;
        id            = createSingleBondEntities(offset_1, col_1, offset_2, col_2, entities, radius / 2.0);
        v_offset      = v_offset * glm::angleAxis((float)(2.0 / 3.0 * glm::pi<float>()), glm::normalize(pos_2 - pos_1));
    }
    return id;
}

//
// Mol file loading and entity creation
//
bool loadMolFile(MolFile &data, std::string_view path) {
    std::ifstream stream(path.data());
    if (!stream) {
        std::cerr << "Error in opening .mol file " << path << ".\n";
        return false;
    }

#ifdef VERBOSE
    std::cout << "Loading .mol file " << path << "\n";
#endif
    return loadMolStream(data, stream);
}

// @todo error checking
bool loadMolStream(MolFile &data, std::istream &stream) {
    std::getline(stream, data.title);

    // @todo read other metadata
    std::getline(stream, data.metadata);
    std::getline(stream, data.comments);

#ifdef VERBOSE
    std::cout << "Title: " << data.title << "\n";
    std::cout << "Comments: " << data.comments << "\n";
    std::cout << "Metadata: " << data.metadata << "\n";
#endif
    int num_atoms, num_bonds;
    stream >> num_atoms >> num_bonds;
    std::string line;
    std::getline(stream, line);    // Skip rest of line @todo

    data.atoms.resize(num_atoms);
    for (int i = 0; i < num_atoms; ++i) {
        if (!stream)
            return false;
        auto &a = data.atoms[i];

        stream >> a.position.x >> a.position.y >> a.position.z;
        stream >> a.symbol;

        int tmp;                       // @todo
        stream >> a.mass_difference >> tmp >> a.hydrogen_count >> tmp >> a.valence;
        std::getline(stream, line);    // Skip rest of line @todo

#ifdef VERBOSE
        std::cout << "ATOM: symbol " << a.symbol << ", position " << a.position << "\n";
#endif
    }

    data.bonds.resize(num_bonds);
    for (int i = 0; i < num_bonds; ++i) {
        if (!stream)
            return false;
        auto &b = data.bonds[i];

        unsigned int type;
        stream >> b.atom_1_index >> b.atom_2_index >> type;
        b.type = (MolBondType)type;
        std::getline(stream, line);    // Skip rest of line @todo

        b.atom_1_index--;
        b.atom_2_index--;
#ifdef VERBOSE
        std::cout << "BOND: number " << i << ", " << b.atom_1_index << " --> " << b.atom_2_index << "\n";
#endif
    }

    return true;
}

void MolFile::clear() {
    atoms.clear();
    bonds.clear();
}

void createEntitiesFromMolFile(Entities &entities, MolFile &data, Camera &camera) {
    static const float relative_cylinder_size = 0.1;
    static const float relative_sphere_size   = 0.26;

    float avg_bond_length = 0.0;
    for (auto &b : data.bonds) {
        auto type = b.type;
        if (b.atom_1_index < 0 || b.atom_1_index > data.atoms.size() ||
            b.atom_2_index < 0 || b.atom_2_index > data.atoms.size())
            continue;
        auto &atom_1 = data.atoms[b.atom_1_index];
        auto &atom_2 = data.atoms[b.atom_2_index];
        avg_bond_length += glm::length(atom_2.position - atom_1.position);
    }
    avg_bond_length /= data.bonds.size();
    float calc_cylinder_r = avg_bond_length * relative_cylinder_size;
    float calc_sphere_r   = avg_bond_length * relative_sphere_size;

    auto center = glm::vec3(0.0);
    entities.instanced_entities.reserve(entities.instanced_entities.size() + data.atoms.size());
    for (auto &a : data.atoms) {
        auto color = getColorFromSymbol(a.symbol);
        createAtomEntity(a.position, color, entities, calc_sphere_r);
        center += a.position;
    }
    center /= data.atoms.size();

    float max_distance = 0.0;
    for (auto &a : data.atoms) {
        auto dist = glm::length(a.position - center);
        if (dist > max_distance) {
            max_distance = dist;
        }
    }

    camera.target_sphere(center, 1.5 * max_distance);

    // Each bond is at least a single bond which is two entities
    entities.instanced_entities.reserve(entities.instanced_entities.size() + data.bonds.size() * 2);
    for (auto &b : data.bonds) {
        auto type = b.type;
        if (b.atom_1_index < 0 || b.atom_1_index > data.atoms.size() ||
            b.atom_2_index < 0 || b.atom_2_index > data.atoms.size())
            continue;
        auto atom_1 = data.atoms[b.atom_1_index];
        auto atom_2 = data.atoms[b.atom_2_index];

        auto color_1 = getColorFromSymbol(atom_1.symbol);
        auto color_2 = getColorFromSymbol(atom_2.symbol);

        switch (type) {
            case MolBondType::SINGLE: {
                createSingleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            case MolBondType::DOUBLE: {
                createDoubleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            case MolBondType::TRIPLE: {
                createTripleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities, calc_cylinder_r);
                break;
            }
            default:
                std::cerr << "Unhandled bond type " << (unsigned int)type << "\n";
                break;
        }
    }
}

//
// Pdb file loading and entity creation
//
long hashBondPair(int id_1, int id_2) {
    if (id_2 < id_1)
        return (long)id_1 ^ ((long)id_2 << 16);
    else
        return (long)id_2 ^ ((long)id_1 << 16);
}

int hashResidueSeqChain(int res_seq, char chain_id) {
    // Since res <= 9999
    return res_seq + 100000 * chain_id;
}

void createResidueId(int seq_num, char res_name[4], char chain_id, char i_code, char result[10]) {
    sprintf(result, "%3s%c%4d%c", res_name, chain_id, seq_num, i_code);
}

glm::vec4 getColorFromResidueName(std::string residue) {
    auto color_1_lu = residue_to_color_lut.find(std::string(residue));
    if (color_1_lu == residue_to_color_lut.end())
        return glm::vec4(unknown_residue_color, 1.0);
    else
        return glm::vec4(color_1_lu->second, 1.0);
}

bool loadPdbDictionaryFile(PdbDictionary &data, std::string_view path) {
    std::ifstream stream(path.data());
    if (!stream) {
        std::cerr << "Error opening .pdb dictionary file " << path << "\n";
        return false;
    }

#ifdef VERBOSE
    std::cout << "Loading .pdb file as dictionary " << path << "\n";
#endif
    return loadPdbDictionaryStream(data, stream);
}

// Dictionary that provides extra information about residue bonds,
// Its format is not the same as a normal pdb file
bool loadPdbDictionaryStream(PdbDictionary &dict, std::istream &stream) {
    char record_name[8];

    std::unordered_map<std::string, PdbDictionaryConnect> *current_residue_connections = nullptr;
    std::string s_line;
    while (stream) {
        std::getline(stream, s_line);
        if (!stream)
            break;
        auto line = s_line.c_str();

        if (strlen(line) >= 6) {
            substrString(line, 0, 6, record_name);

            if (!strcmp(record_name, "RESIDUE")) {
                char res_name[4];
                substrString(line, 10, 12, res_name);

                current_residue_connections = &dict.residues.try_emplace(std::string(res_name)).first->second;

#ifdef VERBOSE
                std::cout << "RESDUE: res_name " << res_name << " \n";
#endif
            } else if (!strcmp(record_name, "CONECT ")) {
                if (current_residue_connections == nullptr)
                    continue;

                char atom_name[5];
                substrString(line, 12, 15, atom_name);

                int num_bonds = 0;
                substrInt(line, 19, 19, &num_bonds);

                char bonded_atom_name[5];

#ifdef VERBOSE
                std::cout << "CONECT: " << atom_name << " num " << num_bonds << "--> ";
#endif

                for (int i = 0; i < num_bonds; i++) {
                    substrString(line, 21 + i * 5, 24 + i * 5, bonded_atom_name);

#ifdef VERBOSE
                    std::cout << bonded_atom_name << ((i == num_bonds - 1) ? "\n" : ", ");
#endif

                    // @todo make this hashing system better
                    std::string atom_name_1 = std::string(atom_name);
                    std::string atom_name_2 = std::string(bonded_atom_name);
                    std::string name_hash;
                    if (atom_name_1.compare(atom_name_2) == -1)    // Compares alphabetically
                        name_hash = atom_name_1 + atom_name_2;
                    else
                        name_hash = atom_name_2 + atom_name_1;
                    auto lu = current_residue_connections->find(name_hash);

                    // Check if a connection exists, if not, create one
                    PdbDictionaryConnect *connect;
                    if (lu == current_residue_connections->end()) {
                        connect              = &current_residue_connections->try_emplace(name_hash).first->second;
                        connect->atom_name_1 = atom_name_1;
                        connect->atom_name_2 = atom_name_2;
                        connect->type        = PdbConnectionType::UNKNOWN;    // 0
                    } else {
                        connect = &lu->second;
                    }
                    connect->type = (PdbConnectionType)glm::min(
                        (unsigned int)connect->type + 1,
                        (unsigned int)PdbConnectionType::TRIPLE_REDUNDANT);    // Fix this C++!
                }
            }
        }
    }

    return true;
}

bool loadPdbFile(PdbFile &data, std::string_view path, const PdbDictionary &dict) {
    std::ifstream stream(path.data());
    if (!stream) {
        std::cerr << "Error in opening .pdb file " << path << ".\n";
        return false;
    }
#ifdef VERBOSE
    std::cout << "Loading .pdb file " << path << "\n";
#endif

    return loadPdbStream(data, stream, dict);
}

bool loadPdbStream(PdbFile &data, std::istream &stream, const PdbDictionary &dict) {
    // @note Assumes line is less than 1024 characters
    char line[1024];
    char record_name[7];

    auto model    = &data.models.try_emplace(0).first->second;
    model->serial = 0;
    // Handle model command being missing
    bool first_model_flag = true;

    std::string s_line;
    while (stream) {
        std::getline(stream, s_line);
        if (!stream)
            break;
        auto line = s_line.c_str();

        if (strlen(line) >= 6) {
            substrString(line, 0, 5, record_name);

            if (!strcmp(record_name, "MODEL ")) {
                int serial;
                sscanf(line, "MODEL %d", &serial);

                if (first_model_flag) {
                    first_model_flag = false;
                } else {
                    if (data.models.find(serial) != data.models.end()) {
#ifdef VERBOSE
                        std::cerr << "Duplicate model serial number, " << serial << "\n";
#endif
                        continue;
                    }
                    model = &data.models.try_emplace(serial).first->second;
                }

                model->serial = serial;
#ifdef VERBOSE
                std::cout << "MODEL " << model->serial << "\n";
#endif
            } else if (!strcmp(record_name, "HETATM") || !strcmp(record_name, "ATOM  ")) {
                int serial;
                substrSscanf(line, 6, 10, " %d", &serial);

                // @note this overwrites duplicate ids
                auto lu0 = model->atoms.find(serial);
#ifdef VERBOSE
                if (lu0 != model->atoms.end())
                    std::cerr << "Collision for serial " << serial << "\n";
#endif

                PdbAtom &atom     = model->atoms[serial];
                atom.serial       = serial;
                atom.is_heterogen = !strcmp(record_name, "HETATM");

                substrString(line, 13, 16, atom.name);

                substrString(line, 17, 19, atom.res_name);
                substrChar(line, 21, &atom.chain_id);
                substrInt(line, 22, 25, &atom.res_seq);

                auto residue_id = hashResidueSeqChain(atom.res_seq, atom.chain_id);
                auto res_lu     = model->residues.find(residue_id);
                // If residue doesn't exist create it
                if (res_lu == model->residues.end()) {
                    auto &residue = model->residues.try_emplace(residue_id).first->second;

                    memcpy(residue.res_name, atom.res_name, 4);
                    residue.chain_id = atom.chain_id;
                    residue.i_code   = atom.i_code;
                    residue.res_seq  = atom.res_seq;

                    // Add new residue to chain
                    auto chain_lu = model->chains.find(atom.chain_id);
                    // If chain doesn't exist create it
                    // @note doesn't account for empty chain_id i.e ' '
                    if (chain_lu == model->chains.end()) {
                        auto &chain    = model->chains.try_emplace(atom.chain_id).first->second;
                        chain.chain_id = atom.chain_id;
                        // @note pointers to values in unordered_map are stable
                        chain.residues.push_back(&residue);

#ifdef VERBOSE
                        std::cout << "CHAIN : " << chain.chain_id << "\n";
#endif
                    } else {
                        auto &chain = chain_lu->second;
                        chain.residues.push_back(&residue);
                    }

                    residue.atom_ids.push_back(atom.serial);
                    residue.atom_name_id[atom.name] = atom.serial;

#ifdef VERBOSE
                    std::cout << "RESDUE: res_name " << residue.res_name << " chain_id " << residue.chain_id
                              << " i_code " << residue.i_code << " res_seq " << residue.res_seq << "\n";
#endif
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

#ifdef VERBOSE
                std::cout << ((atom.is_heterogen) ? "HETATM: " : "ATOM  : ") << " serial " << atom.serial << " name "
                          << atom.name << " symbol " << atom.symbol << " res_name " << atom.res_name << " chain_id "
                          << atom.chain_id << " i_code " << atom.i_code << " res_seq " << atom.res_seq << "\n";
#endif
            } else if (!strcmp(record_name, "CONECT")) {
                int atom_id;
                int connect_ids[4];

                substrInt(line, 6, 10, &atom_id);
                int matches = 0;
                for (; matches < 3; matches++) {
                    if (substrSscanf(line, 11 + 5 * matches, 11 + 5 * matches + 4, " %d", &connect_ids[matches]) != 1)
                        break;
                }

#ifdef VERBOSE
                std::cout << "CONECT: " << atom_id << " --> ";
#endif
                for (int i = 0; i < matches; ++i) {
#ifdef VERBOSE
                    std::cout << connect_ids[i] << ((i == matches - 1) ? "\n" : ", ");
#endif
                    auto hash = hashBondPair(atom_id, connect_ids[i]);
                    auto lu   = model->connections.find(hash);
                    if (lu != model->connections.end()) {
                        auto &b = lu->second;
                        b.type  = PdbConnectionType((unsigned int)b.type + 1);
                    } else {
                        auto b       = &model->connections[hash];
                        b->atom_1_id = atom_id;
                        b->atom_2_id = connect_ids[i];
                        b->type      = PdbConnectionType::SINGLE;
                    }
                }
            } else if (!strcmp(record_name, "HELIX ")) {
                int serial;
                substrInt(line, 7, 9, &serial);
                auto &helix  = model->helices.try_emplace(serial).first->second;
                helix.serial = serial;

                substrString(line, 11, 13, helix.id_code);

                substrString(line, 15, 17, helix.init_res_name);
                substrChar(line, 19, &helix.init_chain_id);
                substrInt(line, 21, 24, &helix.init_seq_num);
                substrChar(line, 25, &helix.init_i_code);

                substrString(line, 27, 29, helix.end_res_name);
                substrChar(line, 31, &helix.end_chain_id);
                substrInt(line, 33, 36, &helix.end_seq_num);
                substrChar(line, 37, &helix.end_i_code);

                substrInt(line, 38, 39, (int *)&helix.type);
                substrString(line, 40, 69, helix.comment);

                substrInt(line, 71, 75, &helix.seq_length);

#ifdef VERBOSE
                std::cout << "HELIX : id_code " << helix.id_code << " init_seq_num " << helix.init_seq_num
                          << " end_seq_num " << helix.end_seq_num << " seq_length " << helix.seq_length << "\n";
#endif
            } else if (!strcmp(record_name, "SHEET ")) {
                char sheet_id[4];
                substrString(line, 11, 13, sheet_id);

                bool is_first = false;
                auto lu       = model->sheets.find(sheet_id);
                PdbSheet *sheet;
                if (lu == model->sheets.end()) {
                    // Create and modify sheet only first time
                    sheet = &model->sheets.try_emplace(sheet_id).first->second;
                    memcpy(sheet->sheet_id, sheet_id, 4);
                    substrInt(line, 14, 15, &sheet->num_strands);
                    is_first = true;
#ifdef VERBOSE
                    std::cout << "SHEET : sheet_id " << sheet->sheet_id << " num_strands " << sheet->num_strands
                              << "\n";
#endif
                } else {
                    sheet = &model->sheets[sheet_id];
                }
                int strand_id;
                substrInt(line, 7, 9, &strand_id);
                auto &strand     = sheet->strands.try_emplace(strand_id).first->second;
                strand.strand_id = strand_id;
                strand.is_first  = is_first;

                substrString(line, 17, 19, strand.init_res_name);
                substrChar(line, 21, &strand.init_chain_id);
                substrInt(line, 22, 25, &strand.init_seq_num);
                substrChar(line, 26, &strand.init_i_code);

                substrString(line, 28, 30, strand.end_res_name);
                substrChar(line, 32, &strand.end_chain_id);
                substrInt(line, 33, 36, &strand.end_seq_num);
                substrChar(line, 37, &strand.end_i_code);

                substrInt(line, 37, 39, &strand.sense);

                if (!is_first) {
                    substrString(line, 41, 44, strand.cur_atom_name);
                    substrString(line, 45, 47, strand.cur_res_name);
                    substrChar(line, 49, &strand.cur_chain_id);
                    substrInt(line, 50, 53, &strand.cur_seq_num);
                    substrChar(line, 54, &strand.cur_i_code);

                    substrString(line, 56, 59, strand.prev_atom_name);
                    substrString(line, 60, 62, strand.prev_res_name);
                    substrChar(line, 64, &strand.prev_chain_id);
                    substrInt(line, 65, 68, &strand.prev_seq_num);
                    substrChar(line, 69, &strand.prev_i_code);
                }
#ifdef VERBOSE
                std::cout << "STRAND: strand_id " << strand.strand_id << "  init_res_name " << strand.init_res_name
                          << "  init_seq_num " << strand.init_seq_num << "  end_res_name " << strand.end_res_name
                          << "  end_seq_num " << strand.end_seq_num << " \n";
#endif
            } else if (!strcmp(record_name, "ENDMDL")) {
#ifdef VERBOSE
                std::cout << "ENDMDL\n";
#endif
            }
        }
    }

    // Processing after loading from file
    for (auto &pair : data.models) {
        auto &m = pair.second;
        // Label residue type based on secondary structures
        for (auto &j : m.helices) {
            auto &helix = j.second;

            if (helix.init_chain_id != helix.end_chain_id)
                fprintf(stderr, "Error: helix chain init and end not equal");
            char chain_id = helix.init_chain_id;

            for (int seq_num = helix.init_seq_num; seq_num <= helix.end_seq_num; seq_num++) {
                auto residue_id = hashResidueSeqChain(seq_num, chain_id);
                auto lu         = m.residues.find(residue_id);
                if (lu == m.residues.end())
                    continue;

                auto &residue = lu->second;
                residue.type  = PdbResidueType::HELIX;
            }
        }
        for (auto &j : m.sheets) {
            auto &sheet = j.second;

            for (auto &k : sheet.strands) {
                auto &strand = k.second;

#ifdef VERBOSE
                if (strand.init_chain_id != strand.end_chain_id)
                    std::cerr << "Strand chain_id init and end not equal";
#endif
                char chain_id = strand.init_chain_id;

                for (int seq_num = strand.init_seq_num; seq_num <= strand.end_seq_num; seq_num++) {
                    auto residue_id = hashResidueSeqChain(seq_num, chain_id);
                    auto lu         = m.residues.find(residue_id);
                    if (lu == m.residues.end())
                        continue;

                    auto &residue = lu->second;
                    residue.type  = PdbResidueType::STRAND;
                }
            }
        }
    }

    addDictionaryToPdb(data, dict);
    return data.models.size();
}

void addDictionaryToPdb(PdbFile &data, const PdbDictionary &dict) {
    // Create connections for unspecified residues from dictionary
    for (auto &pair : data.models) {
        auto &m = pair.second;
        for (auto &rp : m.residues) {
            auto &residue = rp.second;

            auto res_lu = dict.residues.find(std::string(residue.res_name));

            // The residue isn't in our dictionary
            if (res_lu == dict.residues.end()) {
#ifdef VERBOSE
                std::cout << "RESDUE: res_name " << residue.res_name << " chain_id " << residue.chain_id
                          << " i_code " << residue.i_code << " res_seq " << residue.res_seq
                          << " NOT IN PDB DICTIONARY, skipping.\n";
#endif
                continue;
            }

            auto &residue_connections = res_lu->second;
            for (auto &ap : residue_connections) {
                auto atom_1_lu = residue.atom_name_id.find(ap.second.atom_name_1);
                if (atom_1_lu == residue.atom_name_id.end()) {
#ifdef VERBOSE
                    std::cerr << "ATOM  : "
                              << " name " << ap.second.atom_name_1
                              << " NOT IN PDB DICTIONARY'S RESIDUE, skipping\n";
#endif
                    continue;
                }
                auto atom_2_lu = residue.atom_name_id.find(ap.second.atom_name_2);
                if (atom_2_lu == residue.atom_name_id.end()) {
#ifdef VERBOSE
                    std::cerr << "ATOM  : "
                              << " name " << ap.second.atom_name_2
                              << " NOT IN PDB DICTIONARY'S RESIDUE, skipping\n";
#endif
                    continue;
                }
                auto &atom_1_id = atom_1_lu->second;
                auto &atom_2_id = atom_2_lu->second;

                auto hash = hashBondPair(atom_1_id, atom_2_id);
                auto lu   = m.connections.find(hash);
                if (lu == m.connections.end()) {
                    auto &b     = m.connections.try_emplace(hash).first->second;
                    b.atom_1_id = atom_1_id;
                    b.atom_2_id = atom_2_id;
                    b.type      = ap.second.type;
                }
            }
        }
    }
}

// Mesh which blends between different secondary structures in one polypeptide chain
void createPolypeptideChainsMeshes(PdbChain &chain, std::vector<PeptidePlane> &planes) {
    constexpr float r = 0.25;

    // @todo these should probably be dynamic, @note num_splines has to be multiples of 8
    constexpr int num_splines = 16, num_points_per_spline = 12;
    constexpr float arrow_length = 0.7;

    static Profile<num_splines> circle_pf, ellipse_pf, rectangle_pf, ribbon_pf, arrowbase_pf;
    circle_pf.Circle(r);
    ellipse_pf.Ellipse(3 * r, r / 2.f);
    rectangle_pf.Rectangle(6.f * r, 1.5f * r);
    ribbon_pf.Rectangle(5.f * r, 0.5f * r);
    arrowbase_pf.Rectangle(10 * r, 1.5 * r);

    static Shape<4> rectangle_sp, ribbon_sp, arrowbase_sp, arrowtip_sp;
    rectangle_sp.Rectangle(6 * r, 1.5 * r);
    ribbon_sp.Rectangle(5 * r, 0.5 * r);
    arrowbase_sp.Rectangle(10 * r, 1.5 * r);
    arrowtip_sp.Rectangle(0, r);

    // @todo Handle small number of planes
    if (planes.size() < 3)
        return;

    // @todo add less than to template to allow one big buffer
    Spline<num_points_per_spline> spline;
    Spline<2> linear_spline;
    Tube<num_points_per_spline, num_splines, false> tube;
    Tube<2, num_splines, false> linear_tube;
    Tube<num_points_per_spline, 4, true> rectangle_tube;
    Tube<2, 4, true> linear_rectangle_tube;

    Face<num_splines> face;
    Face<4> rectangle_face;

    const int num_planes = planes.size();

    glm::vec3 &prev_n = planes[0].frame.normal;
    for (int i = 0; i < num_planes - 1; i++) {
        auto &p0   = planes[i];
        auto &p1   = planes[i + 1];
        auto &mesh = p0.residue_1->mesh;

        // @note middle peptide or similar might reduce floating point errors
        // need to offset all vertex positions then

        // @todo make parts of mesh different colors
        mesh.num_submeshes = 1;
        mesh.draw_start    = reinterpret_cast<decltype(mesh.draw_start)>(malloc(sizeof(*mesh.draw_start) * mesh.num_submeshes));
        mesh.draw_count    = reinterpret_cast<decltype(mesh.draw_count)>(malloc(sizeof(*mesh.draw_count) * mesh.num_submeshes));
        mesh.draw_start[0] = 0;
        mesh.draw_count[0] = 0;

        Profile<num_splines> *pf0 = &circle_pf, *pf1 = &circle_pf;
        // Flag set when rectangle_sp profiles which aren't transitioning
        const bool do_simple_profile = (p0.residue_1->type == PdbResidueType::HELIX && p1.residue_1->type == PdbResidueType::HELIX) ||
                                       (p0.residue_1->type == PdbResidueType::STRAND && p1.residue_1->type == PdbResidueType::STRAND);
        switch (p0.residue_1->type) {
            case PdbResidueType::HELIX: {
                pf0 = &ribbon_pf;
                break;
            }
            case PdbResidueType::STRAND: {
                pf0 = &rectangle_pf;
                break;
            }
        }
        switch (p1.residue_1->type) {
            case PdbResidueType::HELIX: {
                pf1 = &ribbon_pf;
                break;
            }
            case PdbResidueType::STRAND: {
                pf1 = &rectangle_pf;
                break;
            }
        }

        if (p0.residue_1->type == PdbResidueType::STRAND &&
            p1.residue_1->type == PdbResidueType::COIL) {
            // Normal rectangular spline, but end early to give space for arrow
            spline.Hermite(p0.frame, p1.frame, prev_n, arrow_length);

            rectangle_tube.FromSpline(spline, rectangle_sp, rectangle_sp).add_to_mesh(mesh);

            // Add face on back-side of arrow
            auto bframe = spline.frames[num_points_per_spline - 1];
            rectangle_face.FromShape(arrowbase_sp, bframe).add_to_mesh(mesh, true);

            if (i == num_planes - 2) {
                // Make an arrow tip
                // @hack to make longer arrows at end points
                p1.frame.position += bframe.tangent;

                linear_spline.Hermite(bframe, p1.frame, prev_n);
                linear_rectangle_tube.FromSpline(linear_spline, arrowbase_sp, arrowtip_sp).add_to_mesh(mesh);
            } else {
                // or Attach arrow with next profile
                linear_spline.Hermite(bframe, p1.frame, prev_n);
                linear_tube.FromSpline(linear_spline, arrowbase_pf, *pf1).add_to_mesh(mesh);
            }
        } else {
            // Cap-off beginning of tube
            if (i == 0) {
                // @todo seperate hermite algorithm
                auto &f0 = p0.frame;
                auto &f1 = p1.frame;
                auto f   = f0;

                auto bref = (f0.binormal + f1.binormal) / 2.f;
                f.tangent = glm::normalize(f0.tangent);
                if (glm::length(prev_n) > 0.001) {
                    f.binormal = glm::normalize(glm::cross(f.tangent, prev_n));
                } else {
                    f.binormal = perpendicular_component(bref, f.tangent);
                }
                f.normal = glm::cross(f.binormal, f.tangent);

                face.FromProfile(*pf0, f).add_to_mesh(mesh, true);
            }

            spline.Hermite(p0.frame, p1.frame, prev_n);

            // Rectangle normals are wrong if they are smooth shaded so we need a special case
            if (do_simple_profile) {
                if (p0.residue_1->type == PdbResidueType::HELIX) {
                    rectangle_tube.FromSpline(spline, ribbon_sp, ribbon_sp).add_to_mesh(mesh);
                } else if (p0.residue_1->type == PdbResidueType::STRAND) {
                    rectangle_tube.FromSpline(spline, rectangle_sp, rectangle_sp).add_to_mesh(mesh);
                }
            } else {
                tube.FromSpline(spline, *pf0, *pf1).add_to_mesh(mesh);
            }

            // Cap-off end of tube
            if (i == num_planes - 2) {
                auto &f1   = p1.frame;
                auto f     = f1;
                f.tangent  = glm::normalize(f1.tangent);
                f.normal   = prev_n;
                f.binormal = glm::cross(f.tangent, f.normal);

                face.FromProfile(*pf1, f).add_to_mesh(mesh);
            }
        }

        planes[i].residue_1->mesh_generated = mesh.updateGl();
    }
}

void createPdbModelMeshes(PdbModel &model) {
    for (auto &p : model.chains) {
        auto &chain = p.second;
        if (chain.residues.size() < 2)
            continue;

        std::vector<PeptidePlane> peptide_planes;

        int plane_i = 0;
        for (int i = 0; i < chain.residues.size() - 1; i++) {
            // Determine plane of connection between residues through peptide bond
            auto &r1 = chain.residues[i];
            auto &r2 = chain.residues[i + 1];

            // @note Peptide bond forms between carboxylic group and amino group.
            // PDB should guarantee that residues' label the alpha carbon as "CA"
            // we just use the "O" of the amine to determine the perpendicular plane
            // between adjacent amino acids of protein.
            auto lu_CA1 = r1->atom_name_id.find("CA  ");
            if (lu_CA1 == r1->atom_name_id.end())
                continue;
            auto lu_CA2 = r2->atom_name_id.find("CA  ");
            if (lu_CA2 == r2->atom_name_id.end())
                continue;
            // auto lu_O1  = r1->atom_name_id.find(" O  ");
            // if(lu_O1   == r1->atom_name_id.end()) continue;

            auto CA1 = model.atoms[lu_CA1->second].position;
            auto CA2 = model.atoms[lu_CA2->second].position;
            // auto O1  = model.atoms[lu_O1 ->second].position;

            auto &plane     = peptide_planes.emplace_back();
            plane.residue_1 = r1;
            plane.residue_2 = r2;

            auto &f    = plane.frame;
            f.position = (CA1 + CA2) / 2.f;
            f.tangent  = glm::normalize(CA2 - CA1);
            f.normal   = glm::vec3(0);

            // This a more chemically accurate way of calulating plane frames but is less smooth
            // plane.forward = glm::normalize(CA2 - CA1);
            // plane.normal  = glm::normalize(glm::cross(plane.forward, O1 - CA1));
            //// To ensure perpendicularity calculate right
            // plane.right   = glm::normalize(glm::cross(plane.normal, plane.forward));
        }
        if (peptide_planes.size() == 0)
            continue;

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
            auto &p            = peptide_planes;
            p[i].frame.tangent = (1 - c) * (p[i + 1].frame.position - p[i - 1].frame.position);
        }

        // Uses hermite spline derivatives to determine the normal so adjacent spline's orientations align
        for (uint64_t i = 0; i < (uint64_t)peptide_planes.size() - 1; i++) {
            auto &f0 = peptide_planes[i].frame;
            auto &f1 = peptide_planes[i + 1].frame;

            auto n0 = -6.f * f0.position - 4.f * f0.tangent + 6.f * f1.position - 2.f * f1.tangent;
            n0      = perpendicular_component(n0, f0.tangent);
            f0.normal += n0 / 2.f;

            auto n1 = 6.f * f0.position + 2.f * f0.tangent - 6.f * f1.position + 4.f * f1.tangent;
            n1      = perpendicular_component(n1, f1.tangent);
            f1.normal += n1 / 2.f;
        }

        // Calulcate binormal and ensures adjacent normals aren't flipped
        for (uint64_t i = 0; i < (uint64_t)peptide_planes.size(); i++) {
            auto &f = peptide_planes[i].frame;

            f.normal   = glm::normalize(f.normal);
            f.binormal = glm::normalize(glm::cross(f.tangent, f.normal));
        }
        // This ensures the first and last planes normals are normalised since they only get half a normal
        {
            auto &fb = peptide_planes[0].frame;
            auto &fe = peptide_planes[peptide_planes.size() - 1].frame;

            fb.normal = glm::normalize(glm::cross(fb.tangent, fb.binormal));
            fe.normal = glm::normalize(glm::cross(fe.tangent, fe.binormal));
        }

        // Flip peptides that are >180 from previous
        for (uint64_t i = 1; i < peptide_planes.size(); i++) {
            auto &f0 = peptide_planes[i - 1].frame;
            auto &f  = peptide_planes[i].frame;

            if (glm::dot(f.normal, f0.normal) < 0) {
                f.normal *= -1;
                f.binormal *= -1;
            }
        }

        createPolypeptideChainsMeshes(chain, peptide_planes);
    }
}

// Way of updating entities in-place, definitely needs work @todo make this cleaner
void updateEntityColorsFromPdbModel(Entities &entities, PdbModel &model, PdbDrawSettings &settings) {
    for (const auto &p : model.connections) {
        auto &b    = p.second;
        auto &type = b.type;

        auto atom_1_lu = model.atoms.find(b.atom_1_id);
        if (atom_1_lu == model.atoms.end())
            continue;
        auto atom_1    = atom_1_lu->second;
        auto atom_2_lu = model.atoms.find(b.atom_2_id);
        if (atom_2_lu == model.atoms.end())
            continue;
        auto atom_2 = atom_2_lu->second;

        auto color_1 = getColorFromSymbol(atom_1.symbol);
        auto color_2 = getColorFromSymbol(atom_2.symbol);

        bool is_heterogen_bond = atom_1.is_heterogen && atom_2.is_heterogen;

        if (is_heterogen_bond) {
            if (!settings.draw_hetero_atoms)
                continue;
            color_1.w = settings.hetero_atoms_alpha;
            color_2.w = settings.hetero_atoms_alpha;
        } else {
            if (!settings.draw_residue_atoms)
                continue;
            color_1.w = settings.residue_atoms_alpha;
            color_2.w = settings.residue_atoms_alpha;
        }

        if (b.entity_id != -1) {
            int num_entities = 0;
            switch (type) {
                case PdbConnectionType::SINGLE:
                case PdbConnectionType::SINGLE_REDUNDANT: {
                    num_entities = 2;
                    break;
                }
                case PdbConnectionType::DOUBLE:
                case PdbConnectionType::DOUBLE_REDUNDANT: {
                    num_entities = 4;
                    break;
                }
                case PdbConnectionType::TRIPLE:
                case PdbConnectionType::TRIPLE_REDUNDANT: {
                    num_entities = 6;
                    break;
                }
            }

            bool is_color1 = false;
            for (int id = b.entity_id;
                 id > b.entity_id - num_entities && id >= 0 && id < entities.instanced_entities.size(); --id) {
                if (is_color1)
                    entities.instanced_entities[id].albedo = color_1;
                else
                    entities.instanced_entities[id].albedo = color_2;
                is_color1 = !is_color1;
            }
        }
    }

    for (auto &p : model.atoms) {
        auto &atom = p.second;
        if (atom.entity_id < 0 || atom.entity_id > entities.instanced_entities.size())
            continue;

        auto color = getColorFromSymbol(atom.symbol);
        if (!strcmp(atom.res_name, "HOH")) {
            if (!settings.draw_water_atoms)
                continue;
            color.w = settings.water_atoms_alpha;
        } else if (atom.is_heterogen) {
            if (!settings.draw_hetero_atoms)
                continue;
            color.w = settings.hetero_atoms_alpha;
        } else {
            if (!settings.draw_residue_atoms)
                continue;
            color.w = settings.residue_atoms_alpha;
        }

        entities.instanced_entities[atom.entity_id].albedo = color;
    }

    if (settings.draw_residue_ribbons) {
        for (auto &p : model.chains) {
            auto &chain = p.second;

            auto chain_color = glm::vec4(random_color((unsigned int)chain.chain_id), 1.0);
            for (auto &r : chain.residues) {
                if (r->entity_id < 0 || r->entity_id > entities.mesh_entities.size())
                    continue;

                auto &m_e = entities.mesh_entities[r->entity_id];
                switch (settings.residue_color_mode) {
                    case PdbResidueColorMode::SECONDARY: {
                        switch (r->type) {
                            case PdbResidueType::COIL: {
                                m_e.albedo = coil_color;
                                break;
                            }
                            case PdbResidueType::HELIX: {
                                m_e.albedo = helix_color;
                                break;
                            }
                            case PdbResidueType::STRAND: {
                                m_e.albedo = strand_color;
                                break;
                            }
                        }
                        break;
                    }
                    case PdbResidueColorMode::CHAIN: {
                        m_e.albedo = chain_color;
                        break;
                    }
                    case PdbResidueColorMode::AMINO_ACID: {
                        m_e.albedo = getColorFromResidueName(r->res_name);
                        break;
                    }
                    default: {
#ifdef VERBOSE
                        std::cerr << "Unknown residue color mode" << settings.residue_color_mode << "\n";
#endif
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
    float avg_bond_length                     = 0.0;

    // We need to calculate average bond length before hand to determine size,
    // perhaps this could be done by modifying size matrix
    for (const auto &p : model.connections) {
        auto &b = p.second;

        auto atom_1_lu = model.atoms.find(b.atom_1_id);
        if (atom_1_lu == model.atoms.end())
            continue;
        auto atom_1    = atom_1_lu->second;
        auto atom_2_lu = model.atoms.find(b.atom_2_id);
        if (atom_2_lu == model.atoms.end())
            continue;
        auto atom_2 = atom_2_lu->second;

        avg_bond_length += glm::length(atom_2.position - atom_1.position);
    }
    avg_bond_length /= model.connections.size();
    float calc_cylinder_r = glm::max(avg_bond_length * relative_cylinder_size, 0.1f);
    float calc_sphere_r   = glm::max(avg_bond_length * relative_sphere_size, 0.2f);

    for (auto &p : model.connections) {
        auto &b    = p.second;
        auto &type = b.type;

        auto atom_1_lu = model.atoms.find(b.atom_1_id);
        if (atom_1_lu == model.atoms.end())
            continue;
        auto atom_1    = atom_1_lu->second;
        auto atom_2_lu = model.atoms.find(b.atom_2_id);
        if (atom_2_lu == model.atoms.end())
            continue;
        auto atom_2 = atom_2_lu->second;

        auto color_1 = getColorFromSymbol(atom_1.symbol);
        auto color_2 = getColorFromSymbol(atom_2.symbol);

        bool is_heterogen_bond = atom_1.is_heterogen && atom_2.is_heterogen;

        if (is_heterogen_bond) {
            if (!settings.draw_hetero_atoms)
                continue;
            color_1.w = settings.hetero_atoms_alpha;
            color_2.w = settings.hetero_atoms_alpha;
        } else {
            if (!settings.draw_residue_atoms)
                continue;
            color_1.w = settings.residue_atoms_alpha;
            color_2.w = settings.residue_atoms_alpha;
        }

        switch (type) {
            case PdbConnectionType::SINGLE:
            case PdbConnectionType::SINGLE_REDUNDANT: {
                b.entity_id = createSingleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities,
                                                       calc_cylinder_r);
                break;
            }
            case PdbConnectionType::DOUBLE:
            case PdbConnectionType::DOUBLE_REDUNDANT: {
                b.entity_id = createDoubleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities,
                                                       calc_cylinder_r);
                break;
            }
            case PdbConnectionType::TRIPLE:
            case PdbConnectionType::TRIPLE_REDUNDANT: {
                b.entity_id = createTripleBondEntities(atom_1.position, color_1, atom_2.position, color_2, entities,
                                                       calc_cylinder_r);
                break;
            }
            default:
#ifdef VERBOSE
                std::cerr << "Unknown bond type" << (unsigned int)type << "\n";
#endif
                break;
        }
    }
    auto center = glm::dvec3(0.0);
    for (auto &p : model.atoms) {
        auto &atom = p.second;
        center += glm::dvec3(atom.position) / (double)model.atoms.size();

        auto color = getColorFromSymbol(atom.symbol);
        if (!strcmp(atom.res_name, "HOH")) {
            if (!settings.draw_water_atoms)
                continue;
            color.w = settings.water_atoms_alpha;
        } else if (atom.is_heterogen) {
            if (!settings.draw_hetero_atoms)
                continue;
            color.w = settings.hetero_atoms_alpha;
        } else {
            if (!settings.draw_residue_atoms)
                continue;
            color.w = settings.residue_atoms_alpha;
        }

        atom.entity_id = createAtomEntity(atom.position, color, entities, calc_sphere_r);
    }

    // There is probably a faster algorithm for determining a bounding sphere but for now
    // just loop to find center then loop to find max distance from center
    float max_distance = 0.0;
    for (const auto &p : model.atoms) {
        auto &a   = p.second;
        auto dist = glm::length(glm::dvec3(a.position) - center);
        if (dist > max_distance) {
            max_distance = dist;
        }
    }

    camera.target_sphere(center, max_distance);

    // Create entities for ribbon_sp meshes which have been precalculated
    if (settings.draw_residue_ribbons) {
        for (auto &p : model.chains) {
            auto &chain = p.second;

            auto chain_color = glm::vec4(random_color((unsigned int)chain.chain_id), 1.0);
            for (auto &r : chain.residues) {
                if (!r->mesh_generated)
                    continue;

                r->entity_id = entities.mesh_entities.size();
                auto &m_e    = entities.mesh_entities.emplace_back();

                m_e.mesh = &r->mesh;
                switch (settings.residue_color_mode) {
                    case PdbResidueColorMode::SECONDARY: {
                        switch (r->type) {
                            case PdbResidueType::COIL: {
                                m_e.albedo = coil_color;
                                break;
                            }
                            case PdbResidueType::HELIX: {
                                m_e.albedo = helix_color;
                                break;
                            }
                            case PdbResidueType::STRAND: {
                                m_e.albedo = strand_color;
                                break;
                            }
                            default: {
#ifdef VERBOSE
                                std::cerr << "Unknown residue type when coloring " << (unsigned int)r->type << "\n";
#endif
                                break;
                            }
                        }
                        break;
                    }
                    case PdbResidueColorMode::CHAIN: {
                        m_e.albedo = chain_color;
                        break;
                    }
                    case PdbResidueColorMode::AMINO_ACID: {
                        m_e.albedo = getColorFromResidueName(r->res_name);
                        break;
                    }
                    default: {
#ifdef VERBOSE
                        std::cerr << "Unknown residue color mode" << settings.residue_color_mode << "\n";
#endif
                        break;
                    }
                }
                m_e.albedo.w = settings.residue_ribbons_alpha;
            }
        }
    }
    model.renderable = true;
}

void PdbFile::clear() { models.clear(); }
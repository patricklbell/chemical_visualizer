#ifndef ENTITIES_HPP
#define ENTITIES_HPP
#include <vector>

#include <glm/gtc/quaternion.hpp>

#include "assets.hpp"
#include "glm/detail/type_vec.hpp"

enum EntityType {
    ENTITY = 0,
    MESH_ENTITY = 1,
    INSTANCED_ENTITY = 2,
};

struct Entity {
    EntityType type = EntityType::ENTITY;
    int id;

    Entity(int _id) : id(_id){}
};

struct MeshEntity : Entity {
    Mesh mesh;
    glm::vec3 position      = glm::vec3(0.0);
    glm::quat rotation      = glm::quat(0.0,0.0,0.0,1.0);
    glm::mat3 scale         = glm::mat3(1.0);
    glm::vec3 albedo        = glm::vec3(1.0);

    MeshEntity(int _id) : Entity(_id){
        type = (EntityType)(type | EntityType::MESH_ENTITY);
    }
};

enum class InstanceNames : int {
    NONE = -1,
    SPHERE = 0,
    CYLINDER,
    NUM,
};

struct InstancedEntity : Entity {
    // Index into vector of instanced meshes
    InstanceNames instance_mesh = InstanceNames::NONE;

    glm::vec3 position      = glm::vec3(0.0);
    glm::quat rotation      = glm::quat(0.0,0.0,0.0,1.0);
    glm::mat3 scale         = glm::mat3(1.0);
    glm::vec3 albedo        = glm::vec3(1.0);

    InstancedEntity(int _id) : Entity(_id){
        type = (EntityType)(type | EntityType::INSTANCED_ENTITY);
    }
};


struct Entities {
    std::vector<MeshEntity> mesh_entities;
    std::vector<InstancedEntity> instanced_entities;

    void clear();
    ~Entities();
};

#endif

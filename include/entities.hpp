#ifndef ENTITIES_HPP
#define ENTITIES_HPP
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

#include <mesh.hpp>

struct MeshEntity {
    Mesh* mesh         = nullptr;    // Doesn't own mesh
    glm::vec3 position = glm::vec3(0.0);
    glm::quat rotation = glm::quat(0.0, 0.0, 0.0, 1.0);
    glm::mat3 scale    = glm::mat3(1.0);
    glm::vec4 albedo   = glm::vec4(1.0);
};

struct InstancedEntity {
    enum Type : char {
        SPHERE = 0,
        CYLINDER,
        NUM,
    } type = Type::SPHERE;

    glm::vec3 position = glm::vec3(0.0);
    glm::quat rotation = glm::quat(1.0, 0.0, 0.0, 0.0);
    glm::mat3 scale    = glm::mat3(1.0);
    glm::vec4 albedo   = glm::vec4(1.0);
};

struct Entities {
    std::vector<MeshEntity> mesh_entities;
    std::vector<InstancedEntity> instanced_entities;

    void clear();
};

#endif

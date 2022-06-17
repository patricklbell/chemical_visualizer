#ifndef ENTITIES_HPP
#define ENTITIES_HPP
#include <vector>

#include <glm/gtc/quaternion.hpp>

#include "assets.hpp"
#include "glm/detail/type_vec.hpp"

enum EntityType {
    ENTITY = 0,
    MESH_ENTITY = 1,
};

struct Entity {
    EntityType type = EntityType::ENTITY;
    int id;
    Entity(int _id) : id(_id){}
};

struct MeshEntity : Entity {
    Mesh* mesh = nullptr;
    glm::vec3 position      = glm::vec3(0.0);
    glm::quat rotation      = glm::quat(0.0,0.0,0.0,1.0);
    glm::mat3 scale         = glm::mat3(1.0);
    glm::vec3 albedo        = glm::vec3(1.0);

    MeshEntity(int _id) : Entity(_id){
        type = (EntityType)(type | EntityType::MESH_ENTITY);
    }
};

void freeEntity(Entity *e);
void freeEntities(std::vector<Entity *> &entities);

#endif

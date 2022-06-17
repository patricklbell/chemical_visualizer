#include "entities.hpp"

void freeEntity(Entity *e) {
    switch (e->type) {
        case MESH_ENTITY:
            free((MeshEntity*)e);
            break;
        default:
            free(e);
    }
}

void freeEntities(std::vector<Entity *> &entities) {
    for(auto &e : entities) {
        freeEntity(e);
    }
    entities.clear();
}

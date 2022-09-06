#include "entities.hpp"

void Entities::clear() {
    for(auto &e : mesh_entities) {
        e.mesh.free_resources();
    } 
    mesh_entities.clear();
    instanced_entities.clear();
}

Entities::~Entities() {
    clear();
}

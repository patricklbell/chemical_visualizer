#include "entities.hpp"

void Entities::clear() {
    for(auto &e : mesh_entities) {
        e.mesh.free_resources();
    } 
    mesh_entities.clear();
    for(auto &m : instanced_meshes) {
        m.free_resources();
    }
    instanced_meshes.clear();
    instanced_entities.clear();
}

Entities::~Entities() {
    clear();
}

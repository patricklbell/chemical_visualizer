#ifndef ASSETS_HPP
#define ASSETS_HPP

#include <string>

#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

struct Material {
    glm::fvec3 albedo;
} typedef Material;

struct Mesh {
    unsigned short *indices = NULL;
    int             num_materials = 0;
    int             num_vertices = 0;
    int             num_indices = 0;
    glm::fvec3     *vertices = NULL; 
    glm::fvec3     *normals  = NULL;
    GLuint          indices_vbo;
    GLuint 	        vertices_vbo;
    GLuint 	        normals_vbo;
    GLuint          vao;
    GLenum          draw_mode;
    GLenum          draw_type;
    GLint          *draw_start = NULL;
    GLint          *draw_count = NULL;

    // To get around move/swap etc calling destructor, just call it manually
    // @note could implement move, swap etc...
    void free_resources();
};

void createMeshVao(Mesh &mesh);
bool writeMeshFile(const Mesh &mesh, std::string path);
bool readMeshFile(Mesh &mesh, std::string path);
bool loadMeshWithAssimp(Mesh &mesh, std::string path);

#endif

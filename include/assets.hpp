#ifndef ASSETS_HPP
#define ASSETS_HPP

#include <string>

#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include <assimp/scene.h> 

struct Material {
    glm::fvec3 albedo;
} typedef Material;

struct Mesh {
    unsigned short *indices;
    int             num_materials = 0;
    int             num_vertices = 0;
    int             num_indices = 0;
    glm::fvec3*     vertices; 
    glm::fvec3*     normals;
    GLuint          indices_vbo;
    GLuint 	        vertices_vbo;
    GLuint 	        normals_vbo;
    GLuint          vao;
    GLenum          draw_mode;
    GLenum          draw_type;
    GLint          *draw_start;
    GLint          *draw_count;

    ~Mesh();
};

bool writeMeshFile(const Mesh &mesh, std::string path);
bool readMeshFile(Mesh &mesh, std::string path);
bool loadMeshWithAssimp(Mesh &mesh, std::string path);

#endif

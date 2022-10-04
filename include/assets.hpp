#ifndef ASSETS_HPP
#define ASSETS_HPP

#include <string>

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

enum InstancedMeshType : uint64_t
{
    SPHERE = 0,
    CYLINDER = 1,
    NUM_TYPES,
};

struct Material
{
    glm::fvec3 albedo;
} typedef Material;

struct Mesh
{
    unsigned short *indices = NULL;
    int num_materials = 0;
    int num_vertices = 0;
    int num_indices = 0;
    glm::fvec3 *vertices = NULL;
    glm::fvec3 *normals = NULL;
    GLuint indices_vbo;
    GLuint vertices_vbo;
    GLuint normals_vbo;
    GLuint vao;
    GLenum draw_mode;
    GLenum draw_type;
    GLint *draw_start = NULL;
    GLint *draw_count = NULL;

    ~Mesh();
};

void createMeshVao(Mesh &mesh);
bool writeMeshFile(const Mesh &mesh, std::string path);
bool readMeshFile(Mesh &mesh, std::string path);

#endif

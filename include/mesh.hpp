#ifndef MESH_HPP
#define MESH_HPP

#include <string>

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif

#include <glm/glm.hpp>

struct Material {
    glm::fvec3 albedo;
} typedef Material;

struct Mesh {
    int num_submeshes  = 0;
    int num_vertices   = 0;
    int num_indices    = 0;
    unsigned short *indices = NULL;
    glm::fvec3 *vertices    = NULL;
    glm::fvec3 *normals     = NULL;

    GLuint indices_vbo  = GL_FALSE;
    GLuint vertices_vbo = GL_FALSE;
    GLuint normals_vbo  = GL_FALSE;
    GLuint vao          = GL_FALSE;

    GLenum draw_mode = GL_TRIANGLES;
    GLenum draw_type = GL_UNSIGNED_SHORT;

    GLint *draw_start = NULL;
    GLint *draw_count = NULL;

    // Loads mesh from a .mesh, a custom format which is a binary dump of vertices
    // and normals, for more info see my engine (this is an older version so the
    // write function is not the same, although it is self explanatory)
    bool load(std::string_view path);
    bool updateGl();

    ~Mesh();
};

#endif

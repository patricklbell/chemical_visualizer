#include <mesh.hpp>

#include <string>
#include <stdio.h>
#include <iostream>

#if !defined(EMSCRIPTEN)
#include <GL/glew.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif

#include <glm/glm.hpp>

Mesh::~Mesh() {
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &indices_vbo);
    glDeleteBuffers(1, &vertices_vbo);
    glDeleteBuffers(1, &normals_vbo);

    free(indices);
    free(vertices);
    free(normals);

    free(draw_start);
    free(draw_count);
};

bool Mesh::updateGl() {
    glGenBuffers(1, &vertices_vbo);
    glGenBuffers(1, &indices_vbo);
    glGenBuffers(1, &normals_vbo);

    glGenVertexArrays(1, &vao);
    // bind the vao for writing vbos
    glBindVertexArray(vao);

    // Load the packed vector data into a VBO
    glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
    glBufferData(GL_ARRAY_BUFFER, num_vertices * sizeof(*vertices), &vertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, 0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
    glBufferData(GL_ARRAY_BUFFER, num_vertices * sizeof(*normals), &normals[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, false, 0, 0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices_vbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_indices * sizeof(*indices), &indices[0], GL_STATIC_DRAW);

    glBindVertexArray(0);    // Unbind the VAO

    return true;
}

enum MeshAttributes {
    MESH_ATTRIBUTES_VERTICES = 0,
    MESH_ATTRIBUTES_NORMALS  = 1,
};

const unsigned int MESH_FILE_VERSION = 2;
bool Mesh::load(std::string_view path) {
    FILE *f = fopen(path.data(), "rb");
    if (!f) {
        std::cerr << "Error opening mesh file " << path << "\n";
        return false;
    }

#ifdef VERBOSE
    std::cout << "Loading .mesh file " << path << "\n";
#endif

    unsigned int version;
    fread(&version, sizeof(version), 1, f);
    if (version != MESH_FILE_VERSION) {
        std::cerr << "Invalid mesh file version " << version << " expected " << MESH_FILE_VERSION << "\n";
        return false;
    }

    fread(&num_indices, sizeof(num_indices), 1, f);

    indices = reinterpret_cast<decltype(indices)>(malloc(sizeof(*indices) * num_indices));
    fread(indices, sizeof(*indices), num_indices, f);

    char attributes;
    fread(&attributes, sizeof(attributes), 1, f);
    if (attributes != (MESH_ATTRIBUTES_VERTICES | MESH_ATTRIBUTES_NORMALS)) {
        std::cerr << "Invalid mesh attributes, this loader only support vertices and normals\n";
        return false;
    }

    fread(&num_vertices, sizeof(num_vertices), 1, f);

    vertices = reinterpret_cast<decltype(vertices)>(malloc(sizeof(*vertices) * num_vertices));
    normals  = reinterpret_cast<decltype(normals)>(malloc(sizeof(*normals) * num_vertices));

    fread(vertices, sizeof(*vertices), num_vertices, f);
    fread(normals, sizeof(*normals), num_vertices, f);

    fread(&num_submeshes, sizeof(num_submeshes), 1, f);

    draw_start = reinterpret_cast<decltype(draw_start)>(malloc(sizeof(*draw_start) * num_submeshes));
    draw_count = reinterpret_cast<decltype(draw_count)>(malloc(sizeof(*draw_count) * num_submeshes));
    fread(draw_start, sizeof(*draw_start), num_submeshes, f);
    fread(draw_count, sizeof(*draw_count), num_submeshes, f);

    fclose(f);

    return updateGl();
}

#include <vector>
#include <stdio.h>
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cinttypes>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>

#include "assets.hpp"
#include "shader.hpp"
#include "graphics.hpp"

void Mesh::free_resources(){
    glDeleteVertexArrays(1, &vao);

    free(indices);   
    free(vertices);   
    free(normals);
   
    free(draw_start);   
    free(draw_count);   
};

enum class MeshAttributes : char {
    VERTICES = 0,
    NORMALS  = 1,
    TANGENTS = 2,
    UVS      = 4,
};

const unsigned int MESH_FILE_VERSION = 2;
// For now dont worry about size of types on different platforms
bool writeMeshFile(const Mesh &mesh, std::string path){
    printf("--------------------Save Mesh %s--------------------\n", path.c_str());

    FILE *f;
    f=fopen(path.c_str(), "wb");

    fwrite(&MESH_FILE_VERSION, sizeof(unsigned int), 1, f);

    fwrite(&mesh.num_indices, sizeof(int), 1, f);
    fwrite(mesh.indices, sizeof(unsigned short), mesh.num_indices, f);

    // @todo For now just assume every mesh has every attribute
    char attributes = (char)MeshAttributes::VERTICES | (char)MeshAttributes::NORMALS;
    fwrite(&attributes, sizeof(char), 1, f);

    fwrite(&mesh.num_vertices, sizeof(int), 1, f);

    fwrite(mesh.vertices, sizeof(glm::fvec3), mesh.num_vertices, f);
    fwrite(mesh.normals,  sizeof(glm::fvec3), mesh.num_vertices, f);

    fwrite(&mesh.num_materials, sizeof(int), 1, f);

    // Write material indice ranges
    fwrite(mesh.draw_start, sizeof(GLint), mesh.num_materials, f);
    fwrite(mesh.draw_count, sizeof(GLint), mesh.num_materials, f);

    fclose(f);
    return true;
}
void createMeshVao(Mesh &mesh){
	glGenBuffers(1, &mesh.vertices_vbo);
	glGenBuffers(1, &mesh.indices_vbo);
	glGenBuffers(1, &mesh.normals_vbo);

	glGenVertexArrays(1, &mesh.vao);
	// bind the vao for writing vbos
	glBindVertexArray(mesh.vao);

	// Load the packed vector data into a VBO
	glBindBuffer(GL_ARRAY_BUFFER, mesh.vertices_vbo);
	glBufferData(GL_ARRAY_BUFFER, mesh.num_vertices * sizeof(glm::fvec3), &mesh.vertices[0], GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, 0);
	glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, mesh.normals_vbo);
	glBufferData(GL_ARRAY_BUFFER, mesh.num_vertices * sizeof(glm::fvec3), &mesh.normals[0], GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, false, 0, 0);
	glEnableVertexAttribArray(1);

   	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.indices_vbo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.num_indices * sizeof(unsigned short), &mesh.indices[0], GL_STATIC_DRAW);

	glBindVertexArray(0); //Unbind the VAO
}

bool readMeshFile(Mesh &mesh, std::string path){
    FILE *f;
    f=fopen(path.c_str(), "rb");
    if (!f) {
        fprintf(stderr, "Error in reading mesh file %s.\n", path.c_str());
        return false;
    }

    printf("----------------Loading Mesh File %s----------------\n", path.c_str());

    unsigned int version;
    fread(&version, sizeof(unsigned int), 1, f);
    if(version!=MESH_FILE_VERSION){
        fprintf(stderr, "Invalid mesh file version %u expected %u\n", version, MESH_FILE_VERSION);
        return false;
    }

    fread(&mesh.num_indices, sizeof(int), 1, f);
    mesh.indices = (unsigned short *)malloc(sizeof(unsigned short)*mesh.num_indices);
    fread(mesh.indices, sizeof(unsigned short), mesh.num_indices, f);
    printf("Num of indices %d\n", mesh.num_indices);

    char attributes;
    fread(&attributes, sizeof(char), 1, f);
    
    // @todo For now just assume all attributes
    fread(&mesh.num_vertices, sizeof(int), 1, f);
    printf("Num of vertices %d\n", mesh.num_vertices);

    mesh.vertices = (glm::fvec3*)malloc(sizeof(glm::fvec3)*mesh.num_vertices);
    mesh.normals  = (glm::fvec3*)malloc(sizeof(glm::fvec3)*mesh.num_vertices);
    fread(mesh.vertices, sizeof(glm::fvec3), mesh.num_vertices, f);
    fread(mesh.normals,  sizeof(glm::fvec3), mesh.num_vertices, f);

    // @hardcoded
	mesh.draw_mode = GL_TRIANGLES;
	mesh.draw_type = GL_UNSIGNED_SHORT;

    fread(&mesh.num_materials, sizeof(int), 1, f);
    printf("Num of materials %d\n", mesh.num_materials);

    mesh.draw_start = (GLint*)malloc(sizeof(GLint) * mesh.num_materials);
    mesh.draw_count = (GLint*)malloc(sizeof(GLint) * mesh.num_materials);
    fread(mesh.draw_start, sizeof(GLint), mesh.num_materials, f);
    fread(mesh.draw_count, sizeof(GLint), mesh.num_materials, f);
    
    fclose(f);

    createMeshVao(mesh);
    return true;
}

#ifdef USE_ASSIMP 

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags
#include "assimp/material.h"
#include "assimp/types.h"

const auto ai_import_flags = 
    aiProcess_JoinIdenticalVertices |
    //aiProcess_Triangulate |
    //aiProcess_GenNormals |
    //aiProcess_CalcTangentSpace |
    //aiProcess_RemoveComponent (remove colors) |
    aiProcess_LimitBoneWeights |
    aiProcess_ImproveCacheLocality |
    //aiProcess_RemoveRedundantMaterials |
    //aiProcess_GenUVCoords |
    aiProcess_SortByPType |
    aiProcess_FindDegenerates |
    aiProcess_FindInvalidData |
    //aiProcess_FindInstances |
    aiProcess_ValidateDataStructure |
    aiProcess_OptimizeMeshes |
    //aiProcess_OptimizeGraph |
    aiProcess_Debone;


bool loadMeshWithAssimp(Mesh &mesh, std::string path){
	printf("--------------------Loading Mesh %s With Assimp--------------------\n", path.c_str());

	Assimp::Importer importer;

	const aiScene* scene = importer.ReadFile(path, ai_import_flags);
	if( !scene) {
		fprintf( stderr, importer.GetErrorString());
		getchar();
		return false;
	}

	// Allocate arrays for each mesh 
	mesh.num_materials = scene->mNumMeshes;
	mesh.draw_count = (GLint*)malloc(mesh.num_materials * sizeof(GLint));
	mesh.draw_start = (GLint*)malloc(mesh.num_materials * sizeof(GLint));

	mesh.draw_mode = GL_TRIANGLES;
	mesh.draw_type = GL_UNSIGNED_SHORT;

    int indice_offset = 0;
	for (int i = 0; i < scene->mNumMeshes; ++i) {
		const aiMesh* ai_mesh = scene->mMeshes[i]; 

		mesh.draw_start[i] = indice_offset;
        indice_offset += 3*ai_mesh->mNumFaces;
		mesh.draw_count[i] = 3*ai_mesh->mNumFaces;
        mesh.num_vertices += ai_mesh->mNumVertices;
        mesh.num_indices += 3*ai_mesh->mNumFaces;

        printf("Loading mesh index %d from face %d ---> %d.\n", i, mesh.draw_start[i], mesh.draw_start[i]+mesh.draw_count[i]-1);
	}

    mesh.vertices = (glm::fvec3*)malloc(sizeof(glm::fvec3)*mesh.num_vertices);
    mesh.normals  = (glm::fvec3*)malloc(sizeof(glm::fvec3)*mesh.num_vertices);
    mesh.indices  = (unsigned short*)malloc(sizeof(unsigned short)*mesh.num_indices);
    int vertices_offset = 0, indices_offset = 0;
    for (int j = 0; j < scene->mNumMeshes; ++j) {
		const aiMesh* ai_mesh = scene->mMeshes[j]; 
		for(unsigned int i=0; i<ai_mesh->mNumVertices; i++){
            mesh.vertices[vertices_offset + i] = glm::fvec3(
                ai_mesh->mVertices[i].x,
                ai_mesh->mVertices[i].y,
                ai_mesh->mVertices[i].z
            );
		}
        if(ai_mesh->mNormals != NULL){
            for(unsigned int i=0; i<ai_mesh->mNumVertices; i++){
                mesh.normals[vertices_offset + i] = glm::fvec3(
                    ai_mesh->mNormals[i].x,
                    ai_mesh->mNormals[i].y,
                    ai_mesh->mNormals[i].z
                );
            }
		}
		
		for (unsigned int i=0; i<ai_mesh->mNumFaces; i++){
			// Assumes the model has only triangles.
			mesh.indices[indices_offset + 3*i] = ai_mesh->mFaces[i].mIndices[0];
			mesh.indices[indices_offset + 3*i + 1] = ai_mesh->mFaces[i].mIndices[1];
			mesh.indices[indices_offset + 3*i + 2] = ai_mesh->mFaces[i].mIndices[2];
		}
        vertices_offset += ai_mesh->mNumVertices;
        indices_offset += ai_mesh->mNumFaces*3;
    }
    createMeshVao(mesh);	

	// The "scene" pointer will be deleted automatically by "importer"
	return true;
}

#endif

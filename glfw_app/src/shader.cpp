#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>

#if defined(EMSCRIPTEN)
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include <emscripten/html5.h>
#endif

#if !defined(EMSCRIPTEN)
#include <glad/glad.h>
#else
#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#endif
#include <GLFW/glfw3.h>

#include "shader.hpp"

namespace shader {
    GLuint basic_program;
    struct BasicUniforms basic_uniforms;

    GLuint basic_instanced_program;
    struct BasicUniforms basic_instanced_uniforms;

    std::string glsl_version;
}

using namespace shader;

GLuint loadShader(std::string vertex_fragment_file_path, std::string macro_prepends="") {
	const char *path = vertex_fragment_file_path.c_str();
    #if DEBUG
	printf("Loading shader %s.\n", path);
    #endif
	const char *fragment_macro = "#define COMPILING_FS 1\n";
	const char *vertex_macro   = "#define COMPILING_VS 1\n";
	
	GLuint vertex_shader_id   = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER);

	FILE *fp = fopen(path, "rb");

	if (fp == NULL) {
		fprintf(stderr, "Can't open shader file %s.\n", path);
		return 0;
	}
	fseek(fp, 0L, SEEK_END);
	int num_bytes = ftell(fp);
	
	// @note adds \0 to fread
	rewind(fp); 
	char* shader_code = (char*)malloc((num_bytes + 1) * sizeof(char));
	if(shader_code == NULL)
		return 0;
	fread(shader_code, sizeof(char), num_bytes, fp);
	fclose(fp);
	shader_code[num_bytes] = '\0';
    
	GLint result = GL_FALSE;
	int info_log_length;

    #if DEBUG
	printf("Compiling and linking shader: %s\n", path);
    #endif

    char *vertex_shader_code[] = {(char*)glsl_version.c_str(), (char*)vertex_macro, (char*)macro_prepends.c_str(), shader_code};

	glShaderSource(vertex_shader_id, 4, vertex_shader_code, NULL);
	glCompileShader(vertex_shader_id);

	glGetShaderiv(vertex_shader_id, GL_INFO_LOG_LENGTH, &info_log_length);
	if ( info_log_length > 0 ){
		GLuint fragment_shader_id = glCreateShader(GL_FRAGMENT_SHADER);
		char *vertex_shader_error_message = (char *)malloc(sizeof(char) * (info_log_length+1));
		glGetShaderInfoLog(vertex_shader_id, info_log_length, NULL, vertex_shader_error_message);
		fprintf(stderr, "Error compiling vertex shader %s:\n%s\n", path, vertex_shader_error_message);
		free(vertex_shader_error_message);
    	free(shader_code);
		return GL_FALSE;
	}

	char *fragment_shader_code[] = {(char*)glsl_version.c_str(), (char*)fragment_macro, (char*)macro_prepends.c_str(), shader_code};

	glShaderSource(fragment_shader_id, 4, fragment_shader_code, NULL);
	glCompileShader(fragment_shader_id);

	glGetShaderiv(fragment_shader_id, GL_INFO_LOG_LENGTH, &info_log_length);
	if ( info_log_length > 0 ){
		char *fragment_shader_error_message = (char *)malloc(sizeof(char) * (info_log_length+1));
		glGetShaderInfoLog(fragment_shader_id, info_log_length, NULL, fragment_shader_error_message);
		fprintf(stderr, "Error compiling fragment shader %s:\n%s\n", path, fragment_shader_error_message);
		free(fragment_shader_error_message);
    	free(shader_code);
		return GL_FALSE;
	}


	GLuint program_id = glCreateProgram();
	glAttachShader(program_id, vertex_shader_id);
	glAttachShader(program_id, fragment_shader_id);
	glLinkProgram(program_id);

	glGetProgramiv(program_id, GL_INFO_LOG_LENGTH, &info_log_length);
	if ( info_log_length > 1 ){
		char *program_error_message = (char *)malloc(sizeof(char) * (info_log_length+1));
		glGetProgramInfoLog(program_id, info_log_length, NULL, program_error_message);
		fprintf(stderr, "Error linking shader %s:\n%s\n", path, program_error_message);
		free(program_error_message);
    	free(shader_code);
		return GL_FALSE;
	}

	glDetachShader(program_id, vertex_shader_id);
	glDetachShader(program_id, fragment_shader_id);
	
	glDeleteShader(vertex_shader_id);
	glDeleteShader(fragment_shader_id);

	return program_id;
}

void loadBasicShader(std::string path){
    // Create and compile our GLSL program from the shaders
    auto tmp = basic_program;
    basic_program = loadShader(path);
    if(basic_program == GL_FALSE) {
        basic_program = tmp;
        return;
    }

    // Grab uniforms to modify during rendering
    basic_uniforms.mvp                         = glGetUniformLocation(basic_program, "mvp");
    basic_uniforms.model                       = glGetUniformLocation(basic_program, "model");
    basic_uniforms.sun_color                   = glGetUniformLocation(basic_program, "sun_color");
    basic_uniforms.sun_direction               = glGetUniformLocation(basic_program, "sun_direction");
    basic_uniforms.camera_position             = glGetUniformLocation(basic_program, "camera_position");
    basic_uniforms.albedo                      = glGetUniformLocation(basic_program, "albedo");
}

void loadBasicInstancedShader(std::string path){
    // Create and compile our GLSL program from the shaders
    auto tmp = basic_instanced_program;
    basic_instanced_program = loadShader(path);
    if(basic_instanced_program == GL_FALSE) {
        basic_instanced_program = tmp;
        return;
    }

    // Grab uniforms to modify during rendering
    basic_instanced_uniforms.mvp                         = glGetUniformLocation(basic_instanced_program, "mvp");
    basic_instanced_uniforms.model                       = glGetUniformLocation(basic_instanced_program, "model");
    basic_instanced_uniforms.sun_color                   = glGetUniformLocation(basic_instanced_program, "sun_color");
    basic_instanced_uniforms.sun_direction               = glGetUniformLocation(basic_instanced_program, "sun_direction");
    basic_instanced_uniforms.camera_position             = glGetUniformLocation(basic_instanced_program, "camera_position");
    basic_instanced_uniforms.albedo                      = glGetUniformLocation(basic_instanced_program, "albedo");
}

void loadShader(const std::string &path, ShaderType type) {
    switch(type) {
        case ShaderType::BASIC:
            loadBasicShader(path);
            return;
        case ShaderType::BASIC_INSTANCED:
            loadBasicInstancedShader(path);
            return;
        default:
#if DEBUG
            std::cerr << "Unhandled shader type in loadShader"
#endif
            return;
    }
}

void deleteShaderPrograms(){
    glDeleteProgram(basic_program);
    glDeleteProgram(basic_instanced_program);
}

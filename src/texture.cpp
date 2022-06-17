#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <filesystem>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

GLuint create1x1Texture(const unsigned char color[3], GLint internal_format=GL_RGB){
	GLuint texture_id;
	glGenTextures(1, &texture_id);
	glBindTexture(GL_TEXTURE_2D, texture_id);
	glTexImage2D(GL_TEXTURE_2D, 0, internal_format, 1, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, &color[0]);
	return texture_id;
}
GLuint loadImage(std::string imagepath, GLint internal_format){
	printf("Loading texture: %s\n", imagepath.c_str());

	int x, y, n;
	unsigned char *data = stbi_load(imagepath.c_str(), &x, &y, &n, 3);

	if(data == NULL){
		fprintf(stderr, "Failed to load texture %s.\n", imagepath.c_str());
		return GL_FALSE;
	} 

	GLuint texture_id;
	glGenTextures(1, &texture_id);
	
	glBindTexture(GL_TEXTURE_2D, texture_id);
	glTexImage2D(GL_TEXTURE_2D, 0, internal_format, x, y, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
	stbi_image_free(data);

	// Poor filtering, or ...
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); 

	// ... nice trilinear filtering ...
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glGenerateMipmap(GL_TEXTURE_2D);

	glBindTexture(GL_TEXTURE_2D, 0);

	return texture_id;
}

//GLuint loadCubemap(std::vector<std::string> filenames){
//    GLuint texture_id;
//    glGenTextures(1, &texture_id);
//    glBindTexture(GL_TEXTURE_CUBE_MAP, texture_id);
//
//    for(unsigned int i = 0; i < filenames.size(); i++)
//    {
//        CImg<unsigned char> src(filenames[i].c_str());
//        int w = src.width();
//        int h = src.height();
//        src.permute_axes("cxyz");
//        glTexImage2D(
//            GL_TEXTURE_CUBE_MAP_POSITIVE_X + i,
//            0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, src.data()
//        );
//    }
//    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
//    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
//    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
//	return texture_id;
//}

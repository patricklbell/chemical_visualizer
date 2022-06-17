#ifndef TEXTURE_HPP
#define TEXTURE_HPP

#include <string>
#include <vector>

#include <GLFW/glfw3.h>

GLuint create1x1Texture(const unsigned char color[3], GLint internal_format=GL_RGB);
GLuint loadImage(std::string imagepath, GLint internal_format);
GLuint loadCubemap(std::vector<std::string> filenames);

#endif

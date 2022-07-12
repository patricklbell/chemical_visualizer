#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <string>
#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include "assets.hpp"
#include "controls.hpp"

void scaleMat3(glm::mat3 &mat, glm::vec3 scl);
glm::mat4x4 createModelMatrix(const glm::vec3 &pos, const glm::quat &rot, const glm::mat3x3 &scl);
void screenPosToWorldRay(glm::ivec2 mouse_position, glm::mat4 view, glm::mat4 projection, glm::vec3 &out_origin, glm::vec3 &out_direction);
bool rayIntersectsTriangleCull(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, double &t, double &u, double &v);
bool rayIntersectsTriangle(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, double &t, double &u, double &v);
bool rayIntersectsTriangleTestCull(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction);
bool rayIntersectsTriangleTest(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction);
bool rayIntersectsMesh(Mesh *mesh, const glm::mat4x4 &transform, const Camera &camera, const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, glm::vec3 &collision_point, glm::vec3 &normal);
bool lineIntersectsPlane(const glm::vec3 &plane_origin, const glm::vec3 &plane_normal,const glm::vec3 &line_origin, const glm::vec3 &line_direction, float &t);
float closestDistanceBetweenLines(const glm::vec3 &l1_origin, const glm::vec3 &l1_direction, const glm::vec3 &l2_origin, const glm::vec3 &l2_direction, float &l1_t, float &l2_t);
float closestDistanceBetweenLineCircle(const glm::vec3 &line_origin, const glm::vec3 &line_direction, const glm::vec3 &circle_center, const glm::vec3 &circle_normal, float circle_radius, glm::vec3& point);
glm::vec3 anyPerpendicular(const glm::vec3 &v);
glm::quat quatAlignAxisToDirection(const glm::vec3 &axis, const glm::vec3 &direction);
int substrSscanf(const char *src, int start, int end, const char *format, void *result);
void substrString(const char *src, int start, int end, char *result);
void substrChar(const char *src, int pos, char *result);
void substrInt(const char *src, int start, int end, int *result);
void substrFloat(const char *src, int start, int end, float *result);
glm::vec3 randomColor();

#endif /* ifndef UTILITIES_HPP */

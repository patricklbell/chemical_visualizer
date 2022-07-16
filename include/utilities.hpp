#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <string>
#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include "assets.hpp"
#include "controls.hpp"

void scaleMat3(glm::mat3 &mat, glm::vec3 scl);
glm::mat4x4 createModelMatrix(const glm::vec3 &pos, const glm::quat &rot, const glm::mat3x3 &scl);
inline glm::vec3 toModelSpace(const glm::vec3 &v);
inline glm::vec3 calculateTriangleNormalCCW(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3);
inline void printVector(const glm::vec3 &v);
void screenPosToWorldRay(glm::ivec2 mouse_position, glm::mat4 view, glm::mat4 projection, glm::vec3 &out_origin, glm::vec3 &out_direction);
bool rayIntersectsTriangleCull(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, double &t, double &u, double &v);
bool rayIntersectsTriangle(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, double &t, double &u, double &v);
bool rayIntersectsTriangleTestCull(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction);
bool rayIntersectsTriangleTest(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction);
bool rayIntersectsMesh(const Mesh &mesh, const glm::mat4x4 &transform, const Camera &camera, const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, glm::vec3 &collision_point, glm::vec3 &normal);
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
void createCubicSpline(const glm::vec3 &x0, const glm::vec3 &x1, const glm::vec3 &x2, const glm::vec3 &x3, int n, glm::vec3 *curve);
void createCircularProfile(const int n, glm::vec2 *points, float r);
void projectPointsOnPlane(int num, glm::vec3 p, glm::vec3 u, glm::vec3 v, glm::vec2 *in_points, glm::vec3 *out_points);
void createClosedFacedFromProfile(glm::vec3 center_point, int num_points, glm::vec3 *profile, Mesh &mesh, const bool flipped);
void createClosedSurfaceFromSplines(int num_splines, int num_points_per_spline, glm::vec3 *surface, Mesh &mesh);
glm::vec3 randomColor();

#endif /* ifndef UTILITIES_HPP */

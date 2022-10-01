#ifndef UTILITIES_HPP
#define UTILITIES_HPP
// @todo Split out functions that aren't really utilities, eg. mesh, spline and profile functions

#include <string>
#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include "assets.hpp"
#include "controls.hpp"
#include "loader.hpp"

void printVector(const glm::vec3 &v);
void printVector(const glm::vec2 &v);

// Some helpful math functions for 3D graphics
void scaleMat3(glm::mat3 &mat, glm::vec3 scl);
glm::mat4x4 createModelMatrix(const glm::vec3 &pos, const glm::quat &rot, const glm::mat3x3 &scl);
inline glm::vec3 toModelSpace(const glm::vec3 &v);
inline glm::vec3 calculateTriangleNormalCCW(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3);

// These aren't really used currently but they could be used to select chemical structures
void screenPosToWorldRay(glm::ivec2 mouse_position, glm::mat4 view, glm::mat4 projection, glm::vec3 &out_origin, glm::vec3 &out_direction);
bool rayIntersectsTriangleCull(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, double &t, double &u, double &v);
bool rayIntersectsTriangle(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, double &t, double &u, double &v);
bool rayIntersectsTriangleTestCull(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction);
bool rayIntersectsTriangleTest(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction);
bool rayIntersectsMesh(const Mesh &mesh, const glm::mat4x4 &transform, const Camera &camera, const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, glm::vec3 &collision_point, glm::vec3 &normal);
bool lineIntersectsPlane(const glm::vec3 &plane_origin, const glm::vec3 &plane_normal,const glm::vec3 &line_origin, const glm::vec3 &line_direction, float &t);

float closestDistanceBetweenLines(const glm::vec3 &l1_origin, const glm::vec3 &l1_direction, const glm::vec3 &l2_origin, const glm::vec3 &l2_direction, float &l1_t, float &l2_t);
float closestDistanceBetweenLineCircle(const glm::vec3 &line_origin, const glm::vec3 &line_direction, const glm::vec3 &circle_center, const glm::vec3 &circle_normal, float circle_radius, glm::vec3& point);
glm::vec3 perpendicularComponent(glm::vec3 a, glm::vec3 b);
glm::vec3 anyPerpendicular(const glm::vec3 &v);
glm::quat quatAlignAxisToDirection(const glm::vec3 &axis, const glm::vec3 &direction);

int substrSscanf(const char *src, int start, int end, const char *format, void *result);
void substrString(const char *src, int start, int end, char *result);
void substrChar(const char *src, int pos, char *result);
void substrInt(const char *src, int start, int end, int *result);
void substrFloat(const char *src, int start, int end, float *result);

// 
// Some fairly hacky color generation
//
glm::vec3 randomColor();
glm::vec3 randomSaturatedColor();
glm::vec3 hsvToRgb(const glm::vec3 &hsv);

//
// Spline generation algorithms, take a list of points (and normals) and fill them with points along spline formed between p1 -> p4
// @note Some of these algorithms don't actually require 4 planes, some don't use the normals of the planes
//
void createCubicSplineBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pf2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 &prev_normal);
void createCubicSplineNormalsBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pfn1, glm::vec2 *pf2, glm::vec2 *pfn2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 *normals_tube, glm::vec3 &prev_normal);
// @unused
void createCubicSpline(const glm::vec3 &x0, const glm::vec3 &x1, const glm::vec3 &x2, const glm::vec3 &x3, const int n, glm::vec3 *curve);

void createCubicBezierSplineNormalsBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pfn1, glm::vec2 *pf2, glm::vec2 *pfn2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 *normals_tube, glm::vec3 &prev_normal);

// Hermite Splines is the only complete implementation
void createHermiteSplineBetweenProfiles(const int num_points_per_spline, const int num_splines,  glm::vec2 *pf1, glm::vec2 *pf2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 &prev_normal);
void createHermiteSplineNormalsBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pfn1, glm::vec2 *pf2, glm::vec2 *pfn2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 *normals_tube, glm::vec3 &prev_normal, bool double_normals=false);
// Returns the frames at the end of the partial spline, this could technically replace the non partial form when t=1.0
PeptidePlane createPartialHermiteSplineNormalsBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pfn1, glm::vec2 *pf2, glm::vec2 *pfn2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 *normals_tube, glm::vec3 &prev_normal, float t, bool double_normals=true);

// 
// Create 'profiles' meaning an ordered series of points along shape
//
void createCircularProfile(const int n, glm::vec2 *points, float r);
void createCircularProfileNormals(const int n, glm::vec2 *normals);
void createEllipseProfile(const int n, glm::vec2 *points, float r1, float r2);
void createEllipseProfileNormals(const int n, glm::vec2 *normals, float r1, float r2);
void createRectangleProfile(const int n, glm::vec2 *points, float w, float h);
void createRectangleProfileNormals(const int n, glm::vec2 *normals, float w, float h);
void createRectangleProfileNormalsCorrect(glm::vec2* normals, float w, float h);

void projectPointsOnPlane(const int n, glm::vec3 p, glm::vec3 u, glm::vec3 v, glm::vec2 *in_points, glm::vec3 *out_points);

// Mesh generation functions, take in points/normals and triangulate into a mesh
// If the function doesn't take a normals then it index each triangle seperately resulting in flat shading
void createClosedFacedFromProfile(glm::vec3 center_point, int num_points, glm::vec3 *profile, Mesh &mesh, const bool flipped);
void createClosedSurfaceFromSplines(int num_splines, int num_points_per_spline, glm::vec3 *surface, Mesh &mesh);
void createClosedSurfaceFromSplinesNormals(int num_splines, int num_points_per_spline, glm::vec3 *surface, glm::vec3 *normals, Mesh &mesh);
void createClosedSurfaceFromSplinesDoubledNormals(int num_splines, int num_points_per_spline, glm::vec3* surface, glm::vec3* normals, Mesh& mesh);

#endif /* ifndef UTILITIES_HPP */

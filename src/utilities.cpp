#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <limits>
#include <vector>
#include <cstdint>
#include <map>
#include <unordered_map>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include "glm/detail/func_geometric.hpp"
#include "glm/detail/type_mat.hpp"
#include "glm/detail/type_vec.hpp"
#include "glm/fwd.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include "assets.hpp"
#include "entities.hpp"
#include "globals.hpp"
#include "graphics.hpp"
#include "loader.hpp"
#include "utilities.hpp"
#include "loader.hpp"

void scaleMat3(glm::mat3 &mat, glm::vec3 scl) {
    mat[0][0] *= scl.x;
    mat[1][1] *= scl.y;
    mat[2][2] *= scl.z;
}

glm::mat4x4 createModelMatrix(const glm::vec3 &pos, const glm::quat &rot, const glm::mat3x3 &scl){
    return glm::translate(glm::mat4x4(1.0), pos) * glm::mat4_cast(rot) * glm::mat4x4(scl);
}

inline glm::vec3 toModelSpace(const glm::vec3 &v) {
    return glm::vec3(-v.x, -v.y, v.z);
}

inline glm::vec3 calculateTriangleNormalCCW(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
    return glm::normalize(glm::cross(p2 - p1, p3 - p1));
}

inline void printVector(const glm::vec3 &v) {
    printf("x: %13.10f, y: %13.10f, z: %13.10f\n", v.x, v.y, v.z);
}
inline void printVector(const glm::vec2 &v) {
    printf("x: %13.10f, y: %13.10f\n", v.x, v.y);
}

// https://cadxfem.org/inf/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf
// Cornell university paper describing ray intersection algorithms
static const double epsilon = 0.000001;
bool rayIntersectsTriangleTest(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction){
    glm::vec3 edge1 = vertices[1] - vertices[0];
    glm::vec3 edge2 = vertices[2] - vertices[0];
    glm::vec3 p = glm::cross(ray_direction, edge2);
    float det = glm::dot(edge1, p);
    if(det < epsilon && det > -epsilon) return false;

    float inv_det = 1.0f / det;
    
    auto t = ray_origin - vertices[0];
    float u = glm::dot(t, p) * inv_det;
    if(u < 0.0 || u > 1.0) return false;

    auto q = glm::cross(t, edge1);
    float v = glm::dot(ray_direction, q) * inv_det;
    if(v < 0.0 || u + v > 1.0) return false;

    return true;
}

bool rayIntersectsTriangleTestCull(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction){
    auto edge1 = vertices[1] - vertices[0];
    auto edge2 = vertices[2] - vertices[0];
    auto p = glm::cross(ray_direction, edge2);
    float det = glm::dot(edge1, p);
    if(det < epsilon) return false;

    auto t = ray_origin - vertices[0];
    float u = glm::dot(t, p);
    if(u < 0.0 || u > det) return false;

    auto q = glm::cross(t, edge1);
    float v = glm::dot(ray_direction, q);
    if(v < 0.0 || u + v > det) return false;

    return true;
}

// Returns whether ray intersects and calculates t, u, v where t is the distance 
// from ray origin to the intersection, u and v are the barycentric coordinate
// of the intersection. Tests both sides of face i.e. direction doesn't matter
bool rayIntersectsTriangle(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, double &t, double &u, double &v){
    glm::vec3 edge1 = vertices[1] - vertices[0];
    glm::vec3 edge2 = vertices[2] - vertices[0];
    glm::vec3 p = glm::cross(ray_direction, edge2);
    double det = glm::dot(edge1, p);
    if(det < epsilon && det > -epsilon) return false;

    float inv_det = 1.0f / det;
    
    auto t_vec = ray_origin - vertices[0];
    u = glm::dot(t_vec, p) * inv_det;
    if(u < 0.0 || u > 1.0) return false;

    auto q = glm::cross(t_vec, edge1);
    v = glm::dot(ray_direction, q) * inv_det;
    if(v < 0.0 || u + v > 1.0) return false;

    t = glm::dot(edge2, q) * inv_det;
    return true;
}

bool rayIntersectsMesh(const Mesh &mesh, const glm::mat4x4 &transform, const Camera &camera, const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, glm::vec3 &collision_point, glm::vec3 &normal){
    bool flag = false;
    float min_collision_distance;
    for(int j = 0; j < mesh.num_indices; j+=3){
        const auto p1 = mesh.vertices[mesh.indices[j]];
        const auto p2 = mesh.vertices[mesh.indices[j+1]];
        const auto p3 = mesh.vertices[mesh.indices[j+2]];
        glm::vec3 triangle[3] = {
            glm::vec3(transform * glm::vec4(p1, 1.0)),
            glm::vec3(transform * glm::vec4(p2, 1.0)),
            glm::vec3(transform * glm::vec4(p3, 1.0))
        };
        double u, v, t;
        if(rayIntersectsTriangle(triangle, ray_origin, ray_direction, t, u, v)){
            auto collision_distance = glm::length((ray_origin + ray_direction*(float)t) - camera.position);
            if(!flag || (collision_distance < min_collision_distance)){
                flag = true;
                min_collision_distance = collision_distance;
                collision_point = ray_origin + ray_direction*(float)t;
                // Use uv coordinates and all 3 normals
                normal = mesh.normals[mesh.indices[j]];
            }
        }
    }
    return flag;
}

// Returns whether ray intersects and calculates t, u, v where t is the distance 
// from ray origin to the intersection, u and v are the barycentric coordinate. Culls
// back face so ray must enter front of triangle ?? where CCW winding order is front
bool rayIntersectsTriangleCull(const glm::vec3 vertices[3], const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, double &t, double &u, double &v){
    auto edge1 = vertices[1] - vertices[0];
    auto edge2 = vertices[2] - vertices[0];
    auto p = glm::cross(ray_direction, edge2);
    float det = glm::dot(edge1, p);
    if(det < epsilon) return false;

    auto t_vec = ray_origin - vertices[0];
    u = glm::dot(t_vec, p);
    if(u < 0.0 || u > det) return false;

    auto q = glm::cross(t_vec, edge1);
    v = glm::dot(ray_direction, q);
    if(v < 0.0 || u + v > det) return false;

    t = glm::dot(edge2, q);
    float inv_det = 1.0f / det;
    t *= inv_det;
    u *= inv_det;
    v *= inv_det;
    return true;
}
bool lineIntersectsPlane(const glm::vec3 &plane_origin, const glm::vec3 &plane_normal,const glm::vec3 &line_origin, const glm::vec3 &line_direction, float &t){
    float pn_ld = glm::dot(plane_normal, line_direction);
    // Line parallel to plane
    if(pn_ld < epsilon && pn_ld > -epsilon) return false;

    // p = lo + t*ld
    // pn dot (p - po) = 0
    // pn dot (lo + t*ld - po) = 0
    // t*(pn | ld) = pn | (po - lo)
    // t = pn | (po - lo) / (pn | ld)
    t = glm::dot(plane_normal, plane_origin - line_origin) / (pn_ld);

    return true;
}

float closestDistanceBetweenLineCircle(const glm::vec3 &line_origin, const glm::vec3 &line_direction, const glm::vec3 &circle_center, const glm::vec3 &circle_normal, float circle_radius, glm::vec3& point)
{
    float t;
    // Check if line intersects circle plane, ie not parallel
    if(lineIntersectsPlane(circle_center, circle_normal, line_origin, line_direction, t))
    {
        // get the ray's intersection point on the plane which
        // contains the circle
        const glm::vec3 on_plane = line_origin + t*line_direction;
        // project that point on to the circle's circumference
        point = circle_center + circle_radius * normalize(on_plane - circle_center);
        return glm::length(on_plane - point);
    } else {
        // As line is parallel to circle this doesnt need to be projected?
        point = circle_radius * normalize(-line_direction);

        return glm::length(line_origin - point);
    }
}
float closestDistanceBetweenLines(const glm::vec3 &l1_origin, const glm::vec3 &l1_direction, const glm::vec3 &l2_origin, const glm::vec3 &l2_direction, float &l1_t, float &l2_t)
{
    const glm::vec3 dp = l2_origin - l1_origin;
    const float v12 = dot(l1_direction, l1_direction);
    const float v22 = dot(l2_direction, l2_direction);
    const float v1v2 = dot(l1_direction, l2_direction);

    const float det = v1v2 * v1v2 - v12 * v22;

    if (glm::abs(det) > FLT_MIN)
    {
        const float inv_det = 1.f / det;

        const float dpv1 = dot(dp, l1_direction);
        const float dpv2 = dot(dp, l2_direction);

        l1_t = inv_det * (v22 * dpv1 - v1v2 * dpv2);
        l2_t = inv_det * (v1v2 * dpv1 - v12 * dpv2);

        return glm::length(dp + l2_direction * l2_t - l1_direction * l1_t);
    }
    else
    {
        const glm::vec3 a = glm::cross(dp, l1_direction);
        return std::sqrt(dot(a, a) / v12);
    }
}

// Returns the normalized component of a perpendicular to b
glm::vec3 perpendicularComponent(glm::vec3 a, glm::vec3 b) {
    auto a_dot_b = glm::dot(a, b);
    if(glm::abs(a_dot_b) <= epsilon) {
        printf("epsilon\n");
        return glm::normalize(anyPerpendicular(a));
    } else {
        return glm::normalize(b - a_dot_b*a);
    }
}

void screenPosToWorldRay(
    glm::ivec2 mouse_position,
    glm::mat4 view,                     // Camera position and orientation
    glm::mat4 projection,               // Camera parameters (ratio, field of view, near and far planes)
    glm::vec3& out_origin,              // Ouput : Origin of the ray. Starts at the near plane, so if you want the ray to start at the camera's position instead, ignore this.
    glm::vec3& out_direction            // Ouput : Direction, in world space, of the ray that goes "through" the mouse.
) {

    // The ray Start and End positions, in Normalized Device Coordinates
    glm::vec4 ray_start_NDC(
        ((float)mouse_position.x / (float)window_width - 0.5f) * 2.0f, // [0,1024] -> [-1,1]
        -((float)mouse_position.y / (float)window_height - 0.5f) * 2.0f, // [0, 768] -> [-1,1]
        -1.0, // The near plane maps to Z=-1 in Normalized Device Coordinates
        1.0f
    );
    glm::vec4 ray_end_NDC(
        ((float)mouse_position.x / (float)window_width - 0.5f) * 2.0f,
        -((float)mouse_position.y / (float)window_height - 0.5f) * 2.0f,
        0.0,
        1.0f
    );
    // Inverse of projection and view to obtain NDC to world
    glm::mat4 vp = glm::inverse(projection * view);
    glm::vec4 ray_start_world = vp * ray_start_NDC;
    ray_start_world/=ray_start_world.w;
    glm::vec4 ray_end_world   = vp * ray_end_NDC;
    ray_end_world  /=ray_end_world.w;


    out_direction = glm::vec3(glm::normalize(ray_end_world - ray_start_world));
    out_origin = glm::vec3(ray_start_world);
}

// Returns non zero vector perpendicular to v
glm::vec3 anyPerpendicular(const glm::vec3 &v){
    if(v.z == 0) return glm::vec3(v.y, -v.x, 0);
    return glm::vec3(0, v.z, -v.y);
}

glm::quat quatAlignAxisToDirection(const glm::vec3 &axis, const glm::vec3 &direction){
    auto rot_axis = glm::cross(axis, direction);
    float d = glm::dot(axis, direction);
    // If axis is parallel to direction rot_axis is not unique
    if (d < epsilon && glm::dot(rot_axis, rot_axis) < epsilon){
        return glm::normalize(glm::quat(0, anyPerpendicular(direction)));
    } else {
        float s = std::sqrt(glm::dot(axis, axis) * glm::dot(direction, direction)) + d;
        return glm::normalize(glm::quat(s, rot_axis));
    }
}

int substrSscanf(const char *src, int start, int end, const char *format, void *result) {
    char *substr = (char*)malloc(end - start + 2);
    memcpy(substr, &src[start], end - start + 1);
    substr[end - start + 1] = '\0';

    int matches = sscanf(substr, format, result);

    free(substr);
    return matches;
}

// @note relies on end - start + 2 minimum chars allocated for result
void substrString(const char *src, int start, int end, char *result) {
    memcpy(result, &src[start], end - start + 1);
    result[end - start + 1] = '\0';
}

// @note relies on end - start + 2 minimum chars allocated for result
void substrChar(const char *src, int pos, char *result) {
    (*result) = src[pos];
}

void substrInt(const char *src, int start, int end, int *result) {
    substrSscanf(src, start, end, " %d", result);
}

void substrFloat(const char *src, int start, int end, float *result) {
    substrSscanf(src, start, end, " %f", result);
}

glm::vec3 randomColor() {
    return glm::normalize(glm::vec3(rand(), rand(), rand()));
}

// https://en.wikipedia.org/wiki/HSL_and_HSV#Converting_to_RGB
glm::vec3 hsvToRgb(const glm::vec3 &hsv) {
    auto &h = hsv.r;
    auto &s = hsv.g;
    auto &v = hsv.b;

    float c = s*v;
    int h_i = (int)glm::floor(h*6) % 6;
    float x = c * (1 - glm::abs(h_i % 2 - 1));

    float r,g,b;
    switch (h_i) {
        case 0: case 5: r = c;   break;
        case 1: case 4: r = x;   break;
        case 2: case 3: r = 0.0; break;
    }
    switch (h_i) {
        case 0: case 3: g = x;   break;
        case 1: case 2: g = c;   break;
        case 4: case 5: g = 0.0; break;
    }
    switch (h_i) {
        case 0: case 1: b = 0.0; break;
        case 2: case 5: b = x;   break;
        case 3: case 4: b = c;   break;
    }

    return glm::vec3(r,g,b);
}

glm::vec3 randomSaturatedColor() {
    return hsvToRgb(glm::vec3((double)rand() / (double)RAND_MAX, 0.3, 0.99));
}

// http://www.paulbourke.net/miscellaneous/interpolation
// Cubic interpolation applied independently to each dimension between x1 and x2
void createCubicSpline(const glm::vec3 &x0, const glm::vec3 &x1, const glm::vec3 &x2, const glm::vec3 &x3, const int n, glm::vec3 *curve) {
    // @note maybe use smoother method
    glm::vec3 a0, a1, a2, a3;
    a0 = x3 - x2 - x0 + x1;
    a1 = x0 - x1 - a0;
    a2 = x2 - x0;
    a3 = x1;

    curve[0] = x1;
    for(int i = 1; i < n-1; i++) {
        auto x = glm::vec3((float)i / (n-1));
        // Component-wise multiplication
        curve[i] = a0*x*x*x + a1*x*x + a2*x +a3;
    }
    curve[n-1] = x2;
}

void createCubicSplineNormalsBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pfn1, glm::vec2 *pf2, glm::vec2 *pfn2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 *normals_tube, glm::vec3 &prev_normal) {
    auto &pp1 = p1.position;
    auto &pp2 = p2.position;
    auto &pp3 = p3.position;
    auto &pp4 = p4.position;

    glm::vec3 a0, a1, a2, a3;
    a0 = -0.5f*pp1 + 1.5f*pp2 - 1.5f*pp3 + 0.5f*pp4;
    a1 = pp1 - 2.5f*pp2 + 2.f*pp3 - 0.5f*pp4;
    a2 = -0.5f*pp1 + 0.5f*pp3;
    a3 = pp2;

    glm::vec2 *pf  = (glm::vec2*)malloc(sizeof(glm::vec2)*num_splines);
    glm::vec2 *pfn = (glm::vec2*)malloc(sizeof(glm::vec2)*num_splines);
    for(int i = 0; i < num_points_per_spline; i++) {
        auto t = (float)i / (num_points_per_spline-1);

        auto pp =                a0*t*t*t   + a1*t*t   + a2*t + a3;
        auto pt = glm::normalize(3.f*a0*t*t + 2.f*a1*t + a2);

        glm::vec3 pbn;
        if(glm::length(prev_normal) != 0) pbn = glm::normalize(glm::cross(pt, prev_normal));
        else {
            auto pa =            6.f*a0*t   + 2.f*a1;
            pbn = glm::normalize(glm::cross(pt, pa));
        }
        auto pn =  glm::normalize(glm::cross(pt, pbn));

        if(glm::dot(pn, prev_normal) < 0) {
            pn *= -1.f;
            pbn *= -1.f;
        }

        for(int i = 0; i < num_splines; ++i) {
            pf[i] =  glm::mix(pf1[i],  pf2[i],  t*t*t);
            pfn[i] = glm::normalize(glm::mix(pfn1[i], pfn2[i], t*t*t));
        }

        projectPointsOnPlane(num_splines, pp, pn, pbn, pf,  &spline_tube [num_splines*i]);
        projectPointsOnPlane(num_splines, pp, pn, pbn, pfn, &normals_tube[num_splines*i]);
        prev_normal = pn;
    }
}

void createCubicSplineBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pf2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 &prev_normal) {
    auto &pp1 = p1.position;
    auto &pp2 = p2.position;
    auto &pp3 = p3.position;
    auto &pp4 = p4.position;

    // @note maybe use smoother method
    glm::vec3 a0, a1, a2, a3;
    a0 = -0.5f*pp1 + 1.5f*pp2 - 1.5f*pp3 + 0.5f*pp4;
    a1 = pp1 - 2.5f*pp2 + 2.f*pp3 - 0.5f*pp4;
    a2 = -0.5f*pp1 + 0.5f*pp3;
    a3 = pp2;

    glm::vec2 *pf = (glm::vec2*)malloc(sizeof(glm::vec2)*num_splines);
    for(int i = 0; i < num_points_per_spline; i++) {
        auto t = (float)i / (num_points_per_spline-1);

        auto pp =                a0*t*t*t   + a1*t*t   + a2*t + a3;
        auto pt = glm::normalize(3.f*a0*t*t + 2.f*a1*t + a2);

        glm::vec3 pbn;
        if(prev_normal.length() != 0) pbn = glm::normalize(glm::cross(pt, prev_normal));
        else {
            auto pa =            6.f*a0*t   + 2.f*a1;
            pbn = glm::normalize(glm::cross(pt, pa));
        }
        auto pn =  glm::normalize(glm::cross(pt, pbn));

        if(glm::dot(pn, prev_normal) < 0) {
            pn *= -1.f;
            pbn *= -1.f;
        }

        for(int i = 0; i < num_splines; ++i) {
            pf[i] = glm::mix(pf1[i], pf2[i], t*t*t);
        }

        projectPointsOnPlane(num_splines, pp, pn, pbn, pf,  &spline_tube [num_splines*i]);
        prev_normal = pn;
    }
}

void createHermiteSplineBetweenProfiles(const int num_points_per_spline, const int num_splines,  glm::vec2 *pf1, glm::vec2 *pf2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 &prev_normal) {
    auto &pp0 = p2.position;
    auto &pp1 = p3.position;

    auto &m0 = p2.forward;
    auto &m1 = p3.forward;

    auto &n0 = p2.normal;
    auto &n1 = p3.normal;

    auto &b0 = p2.right;
    auto &b1 = p3.right;
       
    auto bref = (b0 + b1) / 2.f;
    glm::vec2* pf = (glm::vec2*)malloc(sizeof(glm::vec2) * num_splines);
    glm::vec2* pfn = (glm::vec2*)malloc(sizeof(glm::vec2) * num_splines);
    for(int i = 0; i < num_points_per_spline; i++) {
        auto t1 = (float)i / (num_points_per_spline-1);
        auto t2 = t1*t1;
        auto t3 = t1*t2;

        auto p = (2*t3-3*t2+1)*pp0 
               + (t3-2*t2+t1)*m0
               + (-2*t3+3*t2)*pp1
               + (t3-t2)*m1;

        auto tn = (6*t2-6*t1)*pp0 
                + (3*t2-4*t1+1)*m0
                + (-6*t2+6*t1)*pp1
                + (3*t2-2*t1)*m1;

        tn = glm::normalize(tn);
        //auto bn = perpendicularComponent(bref, tn);
        glm::vec3 bn;
        if(glm::length(prev_normal) > 0.001) {
            bn = glm::normalize(glm::cross(tn, prev_normal));
        } else {
            bn = perpendicularComponent(bref, tn);
        }
        auto n = glm::cross(bn, tn);

        for(int i = 0; i < num_splines; ++i) {
            pf[i] =  glm::mix(pf1[i],  pf2[i],  t1);
        }

        projectPointsOnPlane(num_splines, p, bn, n, pf,  &spline_tube [num_splines*i]);
        prev_normal = n;
    }
}

PeptidePlane createPartialHermiteSplineNormalsBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pfn1, glm::vec2 *pf2, glm::vec2 *pfn2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 *normals_tube, glm::vec3 &prev_normal, float t) {
    auto &pp0 = p2.position;
    auto &pp1 = p3.position;

    auto &m0 = p2.forward;
    auto &m1 = p3.forward;

    auto &n0 = p2.normal;
    auto &n1 = p3.normal;

    auto &b0 = p2.right;
    auto &b1 = p3.right;
       
    auto bref = (b0 + b1) / 2.f;
    glm::vec2* pf = (glm::vec2*)malloc(sizeof(glm::vec2) * num_splines);
    glm::vec2* pfn = (glm::vec2*)malloc(sizeof(glm::vec2) * num_splines);
    for(int i = 0; i < num_points_per_spline; i++) {
        auto t1 = ((float)i / (num_points_per_spline-1))*t;
        auto t2 = t1*t1;
        auto t3 = t*t2;

        auto p = (2*t3-3*t2+1)*pp0 
               + (t3-2*t2+t1)*m0
               + (-2*t3+3*t2)*pp1
               + (t3-t2)*m1;

        auto tn = (6*t2-6*t1)*pp0 
                + (3*t2-4*t1+1)*m0
                + (-6*t2+6*t1)*pp1
                + (3*t2-2*t1)*m1;

        tn = glm::normalize(tn);
        //auto bn = perpendicularComponent(bref, tn);
        glm::vec3 bn;
        if(glm::length(prev_normal) > 0.001) {
            bn = glm::normalize(glm::cross(tn, prev_normal));
        } else {
            bn = perpendicularComponent(bref, tn);
        }
        auto n = glm::cross(bn, tn);

        for(int i = 0; i < num_splines; ++i) {
            pf[i] =  glm::mix(pf1[i],  pf2[i],  t1);
            pfn[i] = glm::normalize(glm::mix(pfn1[i], pfn2[i], t1));
        }

        projectPointsOnPlane(num_splines, p, bn, n, pf,  &spline_tube [num_splines*i]);
        projectPointsOnPlane(num_splines, glm::vec3(0.0f), bn, n, pfn, &normals_tube[num_splines*i]);

        if(i == num_points_per_spline - 1) {
            return PeptidePlane{p3.residue_1, p3.residue_2, p, bn, tn, n};
        }
        prev_normal = n;
    }
    return PeptidePlane();
}


void createHermiteSplineNormalsBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pfn1, glm::vec2 *pf2, glm::vec2 *pfn2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 *normals_tube, glm::vec3 &prev_normal) {
    auto &pp0 = p2.position;
    auto &pp1 = p3.position;

    auto &m0 = p2.forward;
    auto &m1 = p3.forward;

    auto &n0 = p2.normal;
    auto &n1 = p3.normal;

    auto &b0 = p2.right;
    auto &b1 = p3.right;
       
    auto bref = (b0 + b1) / 2.f;
    glm::vec2* pf = (glm::vec2*)malloc(sizeof(glm::vec2) * num_splines);
    glm::vec2* pfn = (glm::vec2*)malloc(sizeof(glm::vec2) * num_splines);
    for(int i = 0; i < num_points_per_spline; i++) {
        auto t1 = (float)i / (num_points_per_spline-1);
        auto t2 = t1*t1;
        auto t3 = t1*t2;

        auto p = (2*t3-3*t2+1)*pp0 
               + (t3-2*t2+t1)*m0
               + (-2*t3+3*t2)*pp1
               + (t3-t2)*m1;

        auto tn = (6*t2-6*t1)*pp0 
                + (3*t2-4*t1+1)*m0
                + (-6*t2+6*t1)*pp1
                + (3*t2-2*t1)*m1;

        tn = glm::normalize(tn);
        //auto bn = perpendicularComponent(bref, tn);
        glm::vec3 bn;
        if(glm::length(prev_normal) > 0.001) {
            bn = glm::normalize(glm::cross(tn, prev_normal));
        } else {
            bn = perpendicularComponent(bref, tn);
        }
        auto n = glm::cross(bn, tn);

        for(int i = 0; i < num_splines; ++i) {
            pf[i] =  glm::mix(pf1[i],  pf2[i],  t1);
            pfn[i] = glm::normalize(glm::mix(pfn1[i], pfn2[i], t1));
        }

        //createDebugCartesian(p, -0.4f*tn, 0.4f*n, 0.4f*bn, entities, 0.05);

        projectPointsOnPlane(num_splines, p, bn, n, pf,  &spline_tube [num_splines*i]);
        projectPointsOnPlane(num_splines, glm::vec3(0.0f), bn, n, pfn, &normals_tube[num_splines*i]);
        prev_normal = n;
    }
}

void createCubicBezierSplineNormalsBetweenProfiles(const int num_points_per_spline, const int num_splines, glm::vec2 *pf1, glm::vec2 *pfn1, glm::vec2 *pf2, glm::vec2 *pfn2, const PeptidePlane &p1, const PeptidePlane &p2, const PeptidePlane &p3, const PeptidePlane &p4, glm::vec3 *spline_tube, glm::vec3 *normals_tube, glm::vec3 &prev_normal) {
    auto pp0 = p2.position;
    auto pp1 = p2.position + p2.forward;
    auto pp2 = p3.position - p3.forward;
    auto pp3 = p3.position;

    glm::vec2* pf = (glm::vec2*)malloc(sizeof(glm::vec2) * num_splines);
    glm::vec2* pfn = (glm::vec2*)malloc(sizeof(glm::vec2) * num_splines);
    for(int i = 0; i < num_points_per_spline; i++) {
        auto t = (float)i / (num_points_per_spline-1);
        auto t2 = t*t;
        auto t3 = t*t2;

        auto p  = (1-t)*(1-t)*(1-t)*pp0+3*(1-t)*(1-t)*t*pp1 + 3*(1-t)*t2*pp2 + t3*pp3;
        auto tn = 3*(1-t)*(1-t)*(pp1-pp0) + 6*(1-t)*t*(pp2-pp1) + 3*t2*(pp3-pp2);
        tn = glm::normalize(tn);

        glm::vec3 bn;
        if(glm::length(prev_normal) > 0.001) {
            bn = glm::normalize(glm::cross(tn, prev_normal));
        } else {
            auto a  = 6*(1-t)*(pp2-2.f*pp1+pp0) + 6*t*(pp3-2.f*pp2+pp1);
            bn = glm::normalize(glm::cross(tn, a));
        }

        auto n =  glm::normalize(glm::cross(bn, tn));

        for(int i = 0; i < num_splines; ++i) {
            pf[i] =  glm::mix(pf1[i],  pf2[i],  t3);
            pfn[i] = glm::normalize(glm::mix(pfn1[i], pfn2[i], t3));
        }

        projectPointsOnPlane(num_splines, p, n, bn, pf,  &spline_tube [num_splines*i]);
        projectPointsOnPlane(num_splines, glm::vec3(0.0f), n, bn, pfn, &normals_tube[num_splines*i]);
        prev_normal = n;
    }
}

void createCircularProfile(const int n, glm::vec2 *points, float r=1.0) {
    for(int i = 0; i < n; i++) {
        // t : [0, 2π]
        float t = ((float)i / n) * 2.0 * PI;
        points[i].x = r*glm::cos(t);
        points[i].y = r*glm::sin(t);
    }
}
void createCircularProfileNormals(const int n, glm::vec2 *normals) {
    // Circle's normals are just it's positions
    createCircularProfile(n, normals, 1);
}

void createEllipseProfile(const int n, glm::vec2 *points, float r1=1.0, float r2=1.0) {
    for(int i = 0; i < n; i++) {
        // t : [0, 2π]
        float t = ((float)i / n) * 2.0 * PI;
        points[i].x = r1*glm::cos(t);
        points[i].y = r2*glm::sin(t);
    }
}
void createEllipseProfileNormals(const int n, glm::vec2 *normals, float r1=1.0, float r2=1.0) {
    createEllipseProfile(n, normals, r1, r2);
}

void createRectangleProfile(const int n, glm::vec2 *points, float w=1.0, float h=1.0) {
    auto c = glm::vec2(w/2, h/2);
    auto wc = glm::vec2(w/2, 0);
    auto hc = glm::vec2(0, h/2);

    int m = n / 8;
    for(int i = 0; i < m; i++) {
        float t = ((float)(i+1) / (m+1));
        points[i] = glm::mix(wc, c, t);
    }
    for(int i = 0; i < m; i++) {
        float t = ((float)i / m);
        points[m+i] = glm::mix(c, hc, t);
    }

    // Reflect first quadrant about y axis
    for(int i = 0; i < 2*m; i++) {
        points[2*m + i] = points[2*m - i - 1];
        points[2*m + i].x *= -1;
    }
    // Reflect upper half about x axis
    for(int i = 0; i < 4*m; i++) {
        points[4*m + i] = points[4*m - i - 1];
        points[4*m + i].y *= -1;
    }
}
void createRectangleProfileNormals(const int n, glm::vec2 *normals, float w=1.0, float h=1.0) {
    auto c = glm::vec2(w/2, h/2);
    auto wc = glm::vec2(w/2, 0);
    auto hc = glm::vec2(0, h/2);

    int m = n / 8;
    for(int i = 0; i < m; i++) {
        normals[i] = glm::vec2(1, 0);
    }
    for(int i = 0; i < m; i++) {
        normals[m+i] = glm::vec2(0, 1);
    }

    // Reflect first quadrant about y axis
    for(int i = 0; i < 2*m; i++) {
        normals[2*m + i] = normals[2*m - i - 1];
        normals[2*m + i].x *= -1;
    }
    // Reflect upper half about x axis
    for(int i = 0; i < 4*m; i++) {
        normals[4*m + i] = normals[4*m - i - 1];
        normals[4*m + i].y *= -1;
    }
}

// Transforms 2D points onto 3D plane
void projectPointsOnPlane(const int n, glm::vec3 p, glm::vec3 u, glm::vec3 v, glm::vec2 *in_points, glm::vec3 *out_points) {
    for(int i = 0; i < n; i++){
        auto &in_p = in_points[i];
        out_points[i] = u*in_p.x + v*in_p.y + p;
    }
}

// @note assumes points lie on a plane i.e a profile
void createClosedFacedFromProfile(glm::vec3 center_point, int num_points, glm::vec3 *profile, Mesh &mesh, const bool flipped=false) {
    int vertex_offset = mesh.num_vertices;
    int index_offset = mesh.num_indices;
    
    int &num_triangles = num_points;
    mesh.num_indices += num_triangles*3;
    mesh.indices = (unsigned short *)realloc(mesh.indices, sizeof(*mesh.indices)*mesh.num_indices);
    if(mesh.indices == NULL) fprintf(stderr, "indices Relloc failed\n");

    // All normals are the same so we can use indices for all triangles 
    mesh.num_vertices += num_points + 1;
    mesh.vertices = (glm::vec3 *)realloc(mesh.vertices, sizeof(*mesh.vertices)*mesh.num_vertices);
    if(mesh.vertices == NULL) fprintf(stderr, "vertices Relloc failed\n");
    mesh.normals =  (glm::vec3 *)realloc(mesh.normals,  sizeof(*mesh.normals)*mesh.num_vertices);
    if(mesh.normals == NULL) fprintf(stderr, "normals Relloc failed\n");

    mesh.draw_count[0] += 3*num_triangles;

    glm::vec3 n;
    if(flipped) n = calculateTriangleNormalCCW(center_point, profile[0], profile[1]);
    else        n = calculateTriangleNormalCCW(center_point, profile[1], profile[0]);
    n = toModelSpace(n);

    mesh.vertices[vertex_offset] = toModelSpace(center_point);
    mesh.normals [vertex_offset] = n;

    for(int i = 0; i < num_triangles; i++) {
        mesh.vertices[vertex_offset + i + 1] = toModelSpace(profile[i]);
        mesh.normals [vertex_offset + i + 1] = n;

        // CCW winding order
        mesh.indices [index_offset + 3*i] = vertex_offset;

        if(flipped) {
            mesh.indices [index_offset + 3*i + 1] = vertex_offset + i + 1;
            mesh.indices [index_offset + 3*i + 2] = vertex_offset + (i + 1) % num_triangles + 1;
        } else {
            mesh.indices [index_offset + 3*i + 1] = vertex_offset + (i + 1) % num_triangles + 1;
            mesh.indices [index_offset + 3*i + 2] = vertex_offset + i + 1;
        }
    }
}

// Save reallocs by doing all splines together
// @note order of splines and points must be correct
// @note doesn't recreate vao after modifying mesh
// surface should be a 2d array of num_splines by num_points_per_spline packed into a 1d array
void createClosedSurfaceFromSplines(int num_splines, int num_points_per_spline, glm::vec3 *surface, Mesh &mesh) {
    int vertex_offset = mesh.num_vertices;
    int index_offset = mesh.num_indices;

    int num_quads = (num_points_per_spline - 1)*(num_splines);
    
    int num_triangles = 2*num_quads;
    mesh.num_indices += 3*num_triangles;
    //printf("Realloc indices size %d num %d\n", (int)sizeof(*mesh.indices), mesh.num_indices); // @debug
    mesh.indices = (unsigned short *)realloc(mesh.indices, sizeof(*mesh.indices)*mesh.num_indices);

    // Since each quad is indexed together
    mesh.num_vertices += 4*num_quads;
    //printf("Realloc vertices size %d num %d\n", (int)sizeof(*mesh.vertices), mesh.num_vertices); // @debug
    mesh.vertices = (glm::vec3 *)realloc(mesh.vertices, sizeof(*mesh.vertices)*mesh.num_vertices);
    mesh.normals =  (glm::vec3 *)realloc( mesh.normals,  sizeof(*mesh.normals)*mesh.num_vertices);

    // @note Assumes everything is one material
    mesh.draw_count[0] += 3*num_triangles;

    int quad = 0;
    for(int i = 0; i < num_splines; i++) {
        for(int j = 0; j < num_points_per_spline - 1; j++) {
            // CCW winding order
            auto &p1 = surface[(i+1)%num_splines + (j+1)*num_splines];
            auto &p2 = surface[(i+1)%num_splines + (j  )*num_splines];
            auto &p3 = surface[i                 + (j  )*num_splines];
            auto &p4 = surface[i                 + (j+1)*num_splines];
            auto n = toModelSpace(calculateTriangleNormalCCW(p1, p2, p3));

            mesh.vertices[vertex_offset + 4*quad    ] = toModelSpace(p1);
            mesh.vertices[vertex_offset + 4*quad + 1] = toModelSpace(p2);
            mesh.vertices[vertex_offset + 4*quad + 2] = toModelSpace(p3);
            mesh.vertices[vertex_offset + 4*quad + 3] = toModelSpace(p4);
            mesh.normals [vertex_offset + 4*quad    ] = n;
            mesh.normals [vertex_offset + 4*quad + 1] = n;
            mesh.normals [vertex_offset + 4*quad + 2] = n;
            mesh.normals [vertex_offset + 4*quad + 3] = n;

            mesh.indices [index_offset  + 6*quad    ] = vertex_offset + 4*quad;
            mesh.indices [index_offset  + 6*quad + 1] = vertex_offset + 4*quad + 1;
            mesh.indices [index_offset  + 6*quad + 2] = vertex_offset + 4*quad + 2;

            mesh.indices [index_offset  + 6*quad + 3] = vertex_offset + 4*quad + 2;
            mesh.indices [index_offset  + 6*quad + 4] = vertex_offset + 4*quad + 3;
            mesh.indices [index_offset  + 6*quad + 5] = vertex_offset + 4*quad;

            quad++;
        }
    }
}

void createClosedSurfaceFromSplinesNormals(int num_splines, int num_points_per_spline, glm::vec3 *surface, glm::vec3 *normals, Mesh &mesh) {
    int vertex_offset = mesh.num_vertices;
    int index_offset = mesh.num_indices;

    int num_quads = (num_points_per_spline - 1)*(num_splines);
    
    int num_triangles = 2*num_quads;
    mesh.num_indices += 3*num_triangles;
    //printf("Realloc indices size %d num %d\n", (int)sizeof(*mesh.indices), mesh.num_indices); // @debug
    mesh.indices = (unsigned short *)realloc(mesh.indices, sizeof(*mesh.indices)*mesh.num_indices);

    mesh.num_vertices += num_points_per_spline*num_splines;
    //printf("Realloc vertices size %d num %d\n", (int)sizeof(*mesh.vertices), mesh.num_vertices); // @debug
    mesh.vertices = (glm::vec3 *)realloc(mesh.vertices, sizeof(*mesh.vertices)*mesh.num_vertices);
    mesh.normals =  (glm::vec3 *)realloc( mesh.normals,  sizeof(*mesh.normals)*mesh.num_vertices);

    mesh.draw_count[0] += 3*num_triangles;

    int quad = 0;
    for(int i = 0; i < num_splines; i++) {
        for(int j = 0; j < num_points_per_spline - 1; j++) {
            // CCW winding order
            auto ip1 = (i+1)%num_splines + (j+1)*num_splines;
            auto ip2 = (i+1)%num_splines + (j  )*num_splines;
            auto ip3 = i                 + (j  )*num_splines;
            auto ip4 = i                 + (j+1)*num_splines;

            int point_indice = i*num_points_per_spline + j;
            mesh.vertices[vertex_offset + point_indice] = toModelSpace(surface[point_indice]);
            mesh.normals [vertex_offset + point_indice] = toModelSpace(normals[point_indice]);

            mesh.indices [index_offset  + 6*quad    ] = vertex_offset + ip1;
            mesh.indices [index_offset  + 6*quad + 1] = vertex_offset + ip2;
            mesh.indices [index_offset  + 6*quad + 2] = vertex_offset + ip3;

            mesh.indices [index_offset  + 6*quad + 3] = vertex_offset + ip3;
            mesh.indices [index_offset  + 6*quad + 4] = vertex_offset + ip4;
            mesh.indices [index_offset  + 6*quad + 5] = vertex_offset + ip1;

            quad++;
        }
        int point_indice = i*num_points_per_spline + num_points_per_spline - 1;
        mesh.vertices[vertex_offset + point_indice] = toModelSpace(surface[point_indice]);
        mesh.normals [vertex_offset + point_indice] = toModelSpace(normals[point_indice]);
    }
}

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
#include "utilities.hpp"
#include "texture.hpp"

void scaleMat3(glm::mat3 &mat, glm::vec3 scl) {
    mat[0][0] *= scl.x;
    mat[1][1] *= scl.y;
    mat[2][2] *= scl.z;
}

glm::mat4x4 createModelMatrix(const glm::vec3 &pos, const glm::quat &rot, const glm::mat3x3 &scl){
    return glm::translate(glm::mat4x4(1.0), pos) * glm::mat4_cast(rot) * glm::mat4x4(scl);
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

bool rayIntersectsMesh(Mesh *mesh, const glm::mat4x4 &transform, const Camera &camera, const glm::vec3 &ray_origin, const glm::vec3 &ray_direction, glm::vec3 &collision_point, glm::vec3 &normal){
    bool flag = false;
    float min_collision_distance;
    for(int j = 0; j < mesh->num_indices; j+=3){
        const auto p1 = mesh->vertices[mesh->indices[j]];
        const auto p2 = mesh->vertices[mesh->indices[j+1]];
        const auto p3 = mesh->vertices[mesh->indices[j+2]];
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
                normal = mesh->normals[mesh->indices[j]];
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

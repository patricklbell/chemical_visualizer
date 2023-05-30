#include <utilities/math.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/epsilon.hpp>

// @todo Check maths for better performance
glm::mat4x4 model_matrix(const glm::vec3& pos, const glm::quat& rot, const glm::mat3x3& scl) {
    return glm::translate(glm::mat4x4(1.0), pos) * glm::mat4_cast(rot) * glm::mat4x4(scl);
}

glm::quat rotate_to(const glm::vec3& from, const glm::vec3& to) {
    auto rot_from = glm::cross(from, to);
    float d       = glm::dot(from, to);
    // If from is parallel to to rot_from is not unique
    if (d < EPSILON && glm::dot(rot_from, rot_from) < EPSILON) {
        return glm::normalize(glm::quat(0, any_perpendicular(to)));
    } else {
        float s =
            std::sqrt(glm::dot(from, from) * glm::dot(to, to)) +
            d;
        return glm::normalize(glm::quat(s, rot_from));
    }
}

void scale_mat3(glm::mat3& mat, glm::vec3 scl) {
    mat[0][0] *= scl.x;
    mat[1][1] *= scl.y;
    mat[2][2] *= scl.z;
}

glm::vec3 any_perpendicular(const glm::vec3& v) {
    if (v.z == 0)
        return glm::vec3(v.y, -v.x, 0);
    return glm::vec3(0, v.z, -v.y);
}

glm::vec3 perpendicular_component(const glm::vec3& a, const glm::vec3& b) {
    auto a_dot_b = glm::dot(a, b);
    if (glm::abs(a_dot_b) <= EPSILON) {
        return glm::normalize(any_perpendicular(a));
    } else {
        return glm::normalize(b - a_dot_b * a);
    }
}

void project(glm::vec2* in, glm::vec3* out, int N,
             const glm::vec3& p, const glm::vec3& u, const glm::vec3& v) {
    for (int i = 0; i < N; i++) {
        out[i] = u * in[i].x + v * in[i].y + p;
    }
}

void project(glm::vec2* in, glm::vec3* out, int N, const Frame& f) {
    project(in, out, N, f.position, f.binormal, f.normal);
}
#ifndef UTILITIES_MATH_HPP
#define UTILITIES_MATH_HPP

#include <ostream>
#include <algorithm>

#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

#define EPSILON 0.000001

struct Frame {
    glm::vec3 position;
    glm::vec3 tangent;
    glm::vec3 normal;
    glm::vec3 binormal;
};

glm::mat4x4 model_matrix(const glm::vec3& pos, const glm::quat& rot, const glm::mat3x3& scl);
glm::quat rotate_to(const glm::vec3& from, const glm::vec3& to);
void scale_mat3(glm::mat3& mat, glm::vec3 scl);
glm::vec3 any_perpendicular(const glm::vec3& v);
glm::vec3 perpendicular_component(const glm::vec3& a, const glm::vec3& b);
void project(glm::vec2* in, glm::vec3* out, int N, const glm::vec3& p, const glm::vec3& u, const glm::vec3& v);
void project(glm::vec2* in, glm::vec3* out, int N, const Frame& f);

//
// GLM Templates
//

// Stream operators for some GLM types
template <typename T, glm::precision P>
std::ostream& operator<<(std::ostream& os, const glm::tvec1<T, P>& v) {
    return os << v.x;
}

template <typename T, glm::precision P>
std::ostream& operator<<(std::ostream& os, const glm::tvec2<T, P>& v) {
    return os << v.x << ", " << v.y;
}

template <typename T, glm::precision P>
std::ostream& operator<<(std::ostream& os, const glm::tvec3<T, P>& v) {
    return os << v.x << ", " << v.y << ", " << v.z;
}

template <typename T, glm::precision P>
std::ostream& operator<<(std::ostream& os, const glm::tvec4<T, P>& v) {
    return os << v.x << ", " << v.y << ", " << v.z << ", " << v.w;
}

template <typename T, glm::precision P>
std::ostream& operator<<(std::ostream& os, const glm::tquat<T, P>& v) {
    return os << v.x << ", " << v.y << ", " << v.z << ", " << v.w;
}

template <typename T, glm::precision P>
std::ostream& operator<<(std::ostream& os, const glm::tmat4x4<T, P>& m) {
    return os << "{ \n\t" << m[0][0] << ", " << m[0][1] << ", " << m[0][2] << ", " << m[0][3] << "\n\t" << m[1][0]
              << ", " << m[1][1] << ", " << m[1][2] << ", " << m[1][3] << "\n\t" << m[2][0] << ", " << m[2][1] << ", "
              << m[2][2] << ", " << m[2][3] << "\n\t" << m[3][0] << ", " << m[3][1] << ", " << m[3][2] << ", "
              << m[3][3] << "\n}";
}

template <typename genType>
GLM_FUNC_QUALIFIER genType linearstep(genType edge0, genType edge1, genType x) {
    GLM_STATIC_ASSERT(std::numeric_limits<genType>::is_iec559, "'linearstep' only accept floating-point inputs");

    return glm::clamp((x - edge0) / (edge1 - edge0), genType(0), genType(1));
}
template <typename T, glm::precision P, template <typename, glm::precision> class vecType>
GLM_FUNC_QUALIFIER vecType<T, P> linearstep(T edge0, T edge1, vecType<T, P> const& x) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'linearstep' only accept floating-point inputs");

    return glm::clamp((x - edge0) / (edge1 - edge0), static_cast<T>(0), static_cast<T>(1));
}
template <typename T, glm::precision P, template <typename, glm::precision> class vecType>
GLM_FUNC_QUALIFIER vecType<T, P> linearstep(vecType<T, P> const& edge0, vecType<T, P> const& edge1,
                                            vecType<T, P> const& x) {
    GLM_STATIC_ASSERT(std::numeric_limits<T>::is_iec559, "'linearstep' only accept floating-point inputs");

    return glm::clamp((x - edge0) / (edge1 - edge0), static_cast<T>(0), static_cast<T>(1));
}

#endif    // UTILITIES_MATH_HPP
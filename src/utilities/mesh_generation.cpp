#include <utilities/mesh_generation.hpp>

glm::vec3 to_mesh_space(const glm::vec3& v) {
    return glm::vec3(-v.x, -v.y, v.z);
}

glm::vec3 normal_CCW(const glm::vec3& p1,
                     const glm::vec3& p2,
                     const glm::vec3& p3) {
    return glm::normalize(glm::cross(p2 - p1, p3 - p1));
}
#ifndef CAMERA_CORE_HPP
#define CAMERA_CORE_HPP

#include <glm/glm.hpp>
#include <utilities/math.hpp>

struct Frustrum {
    float near_plane   = 0.05;
    float far_plane    = 150.0;
    float fov_y        = glm::radians(45.0f);
    float aspect_ratio = 1.0;    // width / height
};

struct Camera {
    enum Type {
        PERSPECTIVE = 0,
        ORTHOGRAPHIC,
    };

    Camera(Frustrum _frustrum = Frustrum(), glm::vec3 _position = glm::vec3(3.0), glm::vec3 _target = glm::vec3(0.0), Type _type = Type::PERSPECTIVE);
    void set_frustrum(Frustrum _frustrum = Frustrum());
    void set_position(glm::vec3 _position);
    void set_target(glm::vec3 _target);
    void set_aspect_ratio(float aspect_ratio = 1.0f);
    void set_type(Type _type = Type::PERSPECTIVE);

    void target_sphere(glm::vec3 origin, float radius, glm::vec3 look = glm::vec3(1, 0, 0));

    bool update();    // Returns true if any update to the matrices occured

    Type type;
    Frustrum frustrum;

    glm::vec3 up = glm::vec3(0, 1, 0);
    glm::vec3 position;
    glm::vec3 target;

    glm::vec3 right;
    glm::vec3 forward;
    glm::vec3 head;

    // These are set by the update method
    glm::mat4 view;
    glm::mat4 projection;
    glm::mat4 vp;
    glm::mat4 inv_vp;

    float max_distance = 100.0; // max distance from target, not enforced!

   private:
    bool view_updated       = false;
    bool projection_updated = false;
};

#endif    // CAMERA_CORE_HPP
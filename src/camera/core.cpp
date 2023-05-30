#include <camera/core.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

Camera::Camera(Frustrum _frustrum, glm::vec3 _position, glm::vec3 _target, Type _type) {
    set_frustrum(_frustrum);
    set_position(_position);
    set_target(_target);
    set_type(_type);
    update();
}

void Camera::set_frustrum(Frustrum _frustrum) {
    frustrum           = _frustrum;
    projection_updated = true;
}

void Camera::set_position(glm::vec3 _position) {
    position     = _position;
    view_updated = true;
}

void Camera::set_target(glm::vec3 _target) {
    target       = _target;
    view_updated = true;
}

void Camera::set_aspect_ratio(float aspect_ratio) {
    frustrum.aspect_ratio = aspect_ratio;
    projection_updated    = true;
}

void Camera::set_type(Type _type) {
    type               = _type;
    projection_updated = true;
}

void Camera::target_sphere(glm::vec3 origin, float radius, glm::vec3 look) {
    float fov_x   = glm::atan(glm::tan(frustrum.fov_y / 2.0f) * frustrum.aspect_ratio) * 2.0f;
    float min_fov = glm::min(fov_x, frustrum.fov_y);

    float d = radius / glm::tan(min_fov / 2.0f);

    d = glm::max(d, frustrum.near_plane + radius);

    set_target(origin);
    set_position(target + look * d);

    frustrum.far_plane = glm::distance(position, target) + radius + 10.0f;    // @hardcoded
    max_distance       = frustrum.far_plane - radius;
    projection_updated = true;
}

bool Camera::update() {
    bool update = false;
    if (view_updated) {
        view         = glm::lookAt(position, target, up);
        right        = glm::vec3(glm::transpose(view)[0]);
        forward      = position - target;
        head         = glm::cross(right, forward);
        view_updated = false;
        projection_updated |= type == Type::ORTHOGRAPHIC;
        update = true;
    }

    if (projection_updated) {
        switch (type) {
            case ORTHOGRAPHIC: {
                auto distance             = glm::length(target - position);
                auto ratio_size_per_depth = glm::atan(frustrum.fov_y);

                auto size_y = ratio_size_per_depth * distance;
                auto size_x = ratio_size_per_depth * frustrum.aspect_ratio * distance;
                projection  = glm::ortho(-size_x, size_x, -size_y, size_y,
                                         -frustrum.near_plane, frustrum.far_plane);
                break;
            }
            case PERSPECTIVE: {
                projection = glm::perspective(frustrum.fov_y, frustrum.aspect_ratio, frustrum.near_plane, frustrum.far_plane);
                break;
            }
        }

        projection_updated = false;
        update             = true;
    }
    if (update) {
        vp     = projection * view;
        inv_vp = glm::inverse(vp);
    }

    return update;
}
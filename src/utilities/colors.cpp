#include <utilities/colors.hpp>

// https://en.wikipedia.org/wiki/HSL_and_HSV#Converting_to_RGB
glm::vec3 hsv_to_rgb(const glm::vec3& hsv) {
    auto& h = hsv.r;
    auto& s = hsv.g;
    auto& v = hsv.b;

    float c = s * v;
    int h_i = (int)glm::floor(h * 6) % 6;
    float x = c * (1 - glm::abs(h_i % 2 - 1));

    float r, g, b;
    switch (h_i) {
        case 0:
        case 5:
            r = c;
            break;
        case 1:
        case 4:
            r = x;
            break;
        case 2:
        case 3:
            r = 0.0;
            break;
    }
    switch (h_i) {
        case 0:
        case 3:
            g = x;
            break;
        case 1:
        case 2:
            g = c;
            break;
        case 4:
        case 5:
            g = 0.0;
            break;
    }
    switch (h_i) {
        case 0:
        case 1:
            b = 0.0;
            break;
        case 2:
        case 5:
            b = x;
            break;
        case 3:
        case 4:
            b = c;
            break;
    }

    return glm::vec3(r, g, b);
}

glm::vec3 random_color(unsigned int seed) {
    srand(seed * 87178291199);    // Multiply by prime to offset seeds that are
                                  // close together, eg 'A', 'B' ...
    return hsv_to_rgb(glm::vec3((double)rand() / (double)RAND_MAX, 0.5, 0.99));
}
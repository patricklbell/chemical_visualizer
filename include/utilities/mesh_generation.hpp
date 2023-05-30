#ifndef UTILITIES_MESH_GENERATION
#define UTILITIES_MESH_GENERATION

#ifndef NDEBUG
#include <iostream>
#endif

#include <glm/glm.hpp>
#include <utilities/math.hpp>
#include <mesh.hpp>


glm::vec3 to_mesh_space(const glm::vec3& v);
glm::vec3 normal_CCW(const glm::vec3& p1,
                     const glm::vec3& p2,
                     const glm::vec3& p3);

// A profile is an enclosed, ordered list points and their corresponding normal
template <size_t N>
struct Profile {
    glm::vec2 points[N];
    glm::vec2 normals[N];

    void Circle(float r) { Ellipse(r, r); }

    void Ellipse(float r1, float r2) {
         for (int i = 0; i < N; i++) {
             // t : [0, 2π]
             float t      = ((float)i / N) * 2.0 * glm::pi<float>();
             points[i].x  = r1 * glm::cos(t);
             points[i].y  = r2 * glm::sin(t);
             normals[i].x = r1 * glm::cos(t);
             normals[i].y = r2 * glm::sin(t);
         }
    }

    void Rectangle(float w, float h) {
         auto c  = glm::vec2(w / 2, h / 2);
         auto wc = glm::vec2(w / 2, 0);
         auto hc = glm::vec2(0, h / 2);

         if (N == 4) {
             points[0] = wc + hc;
             points[1] = -wc + hc;
             points[2] = -wc - hc;
             points[3] = wc - hc;
             return;
         }

         int m = N / 8;
         for (int i = 0; i < m; i++) {
             float t   = ((float)(i + 1) / (m + 1));
             points[i] = glm::mix(wc, c, t);
         }
         for (int i = 0; i < m; i++) {
             float t       = ((float)i / m);
             points[m + i] = glm::mix(c, hc, t);
         }

         // Reflect first quadrant about y axis
         for (int i = 0; i < 2 * m; i++) {
             points[2 * m + i] = points[2 * m - i - 1];
             points[2 * m + i].x *= -1;
         }
         // Reflect uper half about x axis
         for (int i = 0; i < 4 * m; i++) {
             points[4 * m + i] = points[4 * m - i - 1];
             points[4 * m + i].y *= -1;
         }

         for (int i = 0; i < m; i++) {
             normals[i] = glm::vec2(1, 0);
         }
         for (int i = 0; i < m; i++) {
             normals[m + i] = glm::vec2(0, 1);
         }

         // Reflect first quadrant about y axis
         for (int i = 0; i < 2 * m; i++) {
             normals[2 * m + i] = normals[2 * m - i - 1];
             normals[2 * m + i].x *= -1;
         }
         // Reflect upper half about x axis
         for (int i = 0; i < 4 * m; i++) {
             normals[4 * m + i] = normals[4 * m - i - 1];
             normals[4 * m + i].y *= -1;
         }
    }
};

// A shape is a more general form of the profile which has
// doubled normals (each point has two normals, one for each connecting line)
template <size_t N>
struct Shape {
    glm::vec2 points[N];
    glm::vec2 normals[2 * N];

    void Rectangle(float w, float h) {
         static_assert(N == 4);

         auto c    = glm::vec2(w / 2, h / 2);
         auto wc   = glm::vec2(w / 2, 0);
         auto hc   = glm::vec2(0, h / 2);
         points[0] = wc + hc;
         points[1] = -wc + hc;
         points[2] = -wc - hc;
         points[3] = wc - hc;

         normals[0] = glm::vec2(1, 0);
         normals[1] = glm::vec2(0, 1);
         normals[2] = glm::vec2(0, 1);
         normals[3] = glm::vec2(-1, 0);

         normals[4] = glm::vec2(-1, 0);
         normals[5] = glm::vec2(0, -1);
         normals[6] = glm::vec2(0, -1);
         normals[7] = glm::vec2(1, 0);
    }
};

template <size_t N>
struct Spline {
    Frame frames[N];

    Spline& Hermite(Frame start, Frame end, glm::vec3& prev_n, float t = 1.0) {
         auto& p0 = start.position;
         auto& p1 = end.position;

         auto& m0 = start.tangent;
         auto& m1 = end.tangent;

         auto& n0 = start.normal;
         auto& n1 = end.normal;

         auto& b0 = start.binormal;
         auto& b1 = end.binormal;

         auto bref = (b0 + b1) / 2.f;
         for (int i = 0; i < N; i++) {
             auto t1 = ((float)i / (N - 1)) * t;
             auto t2 = t1 * t1;
             auto t3 = t1 * t2;

             auto p = (2 * t3 - 3 * t2 + 1) * p0 + (t3 - 2 * t2 + t1) * m0 + (-2 * t3 + 3 * t2) * p1 + (t3 - t2) * m1;

             auto tn = (6 * t2 - 6 * t1) * p0 + (3 * t2 - 4 * t1 + 1) * m0 + (-6 * t2 + 6 * t1) * p1 + (3 * t2 - 2 * t1) * m1;

             tn = glm::normalize(tn);
             glm::vec3 bn;
             if (glm::length(prev_n) > 0.001) {
                 bn = glm::normalize(glm::cross(tn, prev_n));
             } else {
                 bn = perpendicular_component(bref, tn);
             }
             auto n = glm::cross(bn, tn);

            frames[i].position = p;
            frames[i].tangent  = tn;
            frames[i].normal   = n;
            frames[i].binormal = bn;
            prev_n = n;
         }

        return *this;
    }
};

template <size_t N, size_t M, bool AllowDoubled = true>
struct Tube {
    glm::vec3 points[N*M];
    glm::vec3 normals[AllowDoubled ? 2*N*M : N*M];
    bool doubled = AllowDoubled;

    Tube& FromSpline(Spline<N>& spline, Shape<M>& start, Shape<M>& end) {
        static_assert(AllowDoubled);

        Shape<M> blended;

        for (int i = 0; i < N; ++i) {
            auto t      = (float)i / (N - 1);
            auto& frame = spline.frames[i];

            for (int j = 0; j < M; ++j) {
                blended.points[j] = glm::mix(start.points[j], end.points[j], t);

                blended.normals[2 * j]     = glm::normalize(glm::mix(start.normals[2 * j], end.normals[2 * j], t));
                blended.normals[2 * j + 1] = glm::normalize(glm::mix(start.normals[2 * j + 1], end.normals[2 * j + 1], t));
            }

            project(blended.points, &points[M * i], M, frame.position, frame.binormal, frame.normal);
            project(blended.normals, &normals[2 * M * i], 2 * M, glm::vec3(0), frame.binormal, frame.normal);
        }
        doubled = true;

        return *this;
    }

    Tube& FromSpline(Spline<N>& spline, Profile<M>& start, Profile<M>& end) {
        Profile<M> blended;

        for (int i = 0; i < N; ++i) {
            auto t      = (float)i / (N - 1);
            auto& frame = spline.frames[i];

            for (int j = 0; j < M; ++j) {
                blended.points[j]  = glm::mix(start.points[j], end.points[j], t);
                blended.normals[j] = glm::normalize(glm::mix(start.normals[j], end.normals[j], t));
            }

            project(blended.points, &points[M * i], M, frame.position, frame.binormal, frame.normal);
            project(blended.normals, &normals[M * i], M, glm::vec3(0), frame.binormal, frame.normal);
        }
        doubled = false;

        return *this;
    }

    void add_to_mesh(Mesh& mesh) {
        int vertex_offset = mesh.num_vertices;
        int index_offset  = mesh.num_indices;

        int num_quads = (N - 1) * M;

        int num_triangles = 2 * num_quads;
        mesh.num_indices += 3 * num_triangles;
        mesh.indices = reinterpret_cast<decltype(mesh.indices)>(realloc(mesh.indices, sizeof(*mesh.indices) * mesh.num_indices));

        // Double the vertices since adjacent (in profile) aren't shared
        mesh.num_vertices += (doubled ? 2 : 1) * N * M;
        mesh.vertices = reinterpret_cast<decltype(mesh.vertices)>(realloc(mesh.vertices, sizeof(*mesh.vertices) * mesh.num_vertices));
        mesh.normals  = reinterpret_cast<decltype(mesh.normals)>(realloc(mesh.normals, sizeof(*mesh.normals) * mesh.num_vertices));

        if (mesh.indices == NULL || mesh.vertices == NULL || mesh.normals == NULL) {
#ifndef NDEBUG
            std::cerr << "Reallocs failed in mesh creation, you might have run "
                         "out of memory. Skipping this mesh\n";
#endif
            return;
        }

        mesh.draw_count[0] += 3 * num_triangles;

        int quad = 0;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N - 1; j++) {
                // CCW winding order
                int ip1, ip2, ip3, ip4;
                if (doubled) {
                    ip1 = (2 * (i + 1)) % (2 * M) + 2 * (j + 1) * M;
                    ip2 = (2 * (i + 1)) % (2 * M) + 2 * (j)*M;
                    ip3 = (2 * i + 1) % (2 * M) + 2 * (j)*M;
                    ip4 = (2 * i + 1) % (2 * M) + 2 * (j + 1) * M;
                } else {
                    ip1 = (i + 1) % M + (j + 1) * M;
                    ip2 = (i + 1) % M + (j)*M;
                    ip3 = i + (j)*M;
                    ip4 = i + (j + 1) * M;
                }

                int position_indice = i * N + j;
                int point_indice    = (doubled ? 2 : 1) * position_indice;

                mesh.vertices[vertex_offset + point_indice] = to_mesh_space(points[position_indice]);

                mesh.normals[vertex_offset + point_indice] = to_mesh_space(normals[point_indice]);

                if (doubled) {
                    mesh.vertices[vertex_offset + point_indice + 1] = to_mesh_space(points[position_indice]);
                    mesh.normals[vertex_offset + point_indice + 1]  = to_mesh_space(normals[point_indice + 1]);
                }

                mesh.indices[index_offset + 6 * quad]     = vertex_offset + ip1;
                mesh.indices[index_offset + 6 * quad + 1] = vertex_offset + ip2;
                mesh.indices[index_offset + 6 * quad + 2] = vertex_offset + ip3;

                mesh.indices[index_offset + 6 * quad + 3] = vertex_offset + ip3;
                mesh.indices[index_offset + 6 * quad + 4] = vertex_offset + ip4;
                mesh.indices[index_offset + 6 * quad + 5] = vertex_offset + ip1;

                quad++;
            }

            int position_indice = i * N + N - 1;
            int point_indice    = (doubled ? 2 : 1) * position_indice;

            mesh.vertices[vertex_offset + point_indice] = to_mesh_space(points[position_indice]);

            mesh.normals[vertex_offset + point_indice] = to_mesh_space(normals[point_indice]);

            if (doubled) {
                mesh.vertices[vertex_offset + point_indice + 1] =
                    to_mesh_space(points[position_indice]);
                mesh.normals[vertex_offset + point_indice + 1] =
                    to_mesh_space(normals[point_indice + 1]);
            }
        }
    }
};

template <size_t N>
struct Face {
    glm::vec3 points[N];
    glm::vec3 center = glm::vec3(0);

    Face& FromShape(Shape<N>& sp, Frame& f) {
        project(sp.points, points, N, f);
        center = f.position;
        return *this;
    }
    Face& FromProfile(Profile<N>& pf, Frame& f) {
        project(pf.points, points, N, f);
        center = f.position;
        return *this;
    }

    void add_to_mesh(Mesh& mesh, const bool flipped = false) {
        int vertex_offset = mesh.num_vertices;
        int index_offset  = mesh.num_indices;

        int num_triangles = N;
        mesh.num_indices += num_triangles * 3;
        mesh.indices = reinterpret_cast<decltype(mesh.indices)>(realloc(mesh.indices, sizeof(*mesh.indices) * mesh.num_indices));

        // All normals are the same so we can use indices for all triangles
        mesh.num_vertices += N + 1;
        mesh.vertices = reinterpret_cast<decltype(mesh.vertices)>(realloc(mesh.vertices, sizeof(*mesh.vertices) * mesh.num_vertices));
        mesh.normals  = reinterpret_cast<decltype(mesh.normals)>(realloc(mesh.normals, sizeof(*mesh.normals) * mesh.num_vertices));

        if (mesh.indices == NULL || mesh.vertices == NULL || mesh.normals == NULL) {
#ifndef NDEBUG
            std::cerr << "Reallocs failed in mesh creation, you might have run "
                         "out of memory. Skipping this mesh\n";
#endif
            return;
        }

        mesh.draw_count[0] += 3 * num_triangles;

        glm::vec3 n;
        if (flipped)
            n = normal_CCW(center, points[0], points[1]);
        else
            n = normal_CCW(center, points[1], points[0]);
        n = to_mesh_space(n);

        mesh.vertices[vertex_offset] = to_mesh_space(center);
        mesh.normals[vertex_offset]  = n;

        for (int i = 0; i < num_triangles; i++) {
            mesh.vertices[vertex_offset + i + 1] = to_mesh_space(points[i]);
            mesh.normals[vertex_offset + i + 1]  = n;

            // CCW winding order
            mesh.indices[index_offset + 3 * i] = vertex_offset;

            if (flipped) {
                mesh.indices[index_offset + 3 * i + 1] = vertex_offset + i + 1;
                mesh.indices[index_offset + 3 * i + 2] =
                    vertex_offset + (i + 1) % num_triangles + 1;
            } else {
                mesh.indices[index_offset + 3 * i + 1] =
                    vertex_offset + (i + 1) % num_triangles + 1;
                mesh.indices[index_offset + 3 * i + 2] = vertex_offset + i + 1;
            }
        }
    }
};

#endif    // UTILITIES_MESH_GENERATION
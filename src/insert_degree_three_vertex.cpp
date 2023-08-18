#include "insert_degree_three_vertex.h"

#include <array>

#include "next.h"
#include "twin.h"
#include "glue_face_sides.h"
#include "barycentric_dot_product.h"
#include "is_same_face_side.h"
#include "vertex_one_ring_face_sides.h"
#include "opposite_corner_angle.h"
#include "get_smallest_angular_coordinate.h"
#include "is_boundary_vertex.h"
#include "pi.h"

// borrowed from remove_degree_three_vertex.cpp
static void relabel(
                    const Eigen::Matrix<int, 3, 2> & boundary_fs,
                    const int & f_new,
                    Eigen::Matrix<int, 3, 2> & boundary_fs_twin)
{
    using namespace std;
    using namespace Eigen;

    Vector2i bfs0 = boundary_fs.row(0);
    Vector2i bfs1 = boundary_fs.row(1);
    Vector2i bfs2 = boundary_fs.row(2);
    for (int iS=0; iS<3; iS++) {
        if (is_same_face_side(bfs0, boundary_fs_twin.row(iS))){
            // cout << "relabel: "  << boundary_fs_twin.row(iS).transpose() << f_new << "  0" << endl;
            boundary_fs_twin.row(iS) << f_new, 0;
            continue;
        }
        else if (is_same_face_side(bfs1, boundary_fs_twin.row(iS))){
            // cout << "relabel: "  << boundary_fs_twin.row(iS).transpose() << f_new << "  1" << endl;
            boundary_fs_twin.row(iS) << f_new, 1;
            continue;
        }
        else if (is_same_face_side(bfs2, boundary_fs_twin.row(iS))){
            // cout << "relabel: "  << boundary_fs_twin.row(iS).transpose() << f_new << "  2" << endl;
            boundary_fs_twin.row(iS) << f_new, 2;
            continue;
        }
    }
}

int insert_degree_three_vertex(
                                int f,
                                const Eigen::Vector3d & b,
                                Eigen::MatrixXi & F,
                                Eigen::MatrixXi & G,
                                Eigen::MatrixXd & l,
                                Eigen::MatrixXd & A,
                                Eigen::MatrixXi & v2fs,
                                Eigen::MatrixXd & BC,
                                std::vector<std::vector<int>> & F2V)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    if (f == GHOST_INDEX) {
        throw std::invalid_argument("tried to insert a vertex into a nonexistent face");
    }

    auto inside_triangle = [](const Vector3d& v) -> bool {
        return v(0) >= 0 && v(1) >= 0 && v(2) >= 0;
    };

    if (!inside_triangle(b)) {
        throw std::invalid_argument("vertex insertion point is outside of triangle");
    }

    auto vertex_angle_sum = [&](int iV) -> double {
        vector<Vector2i> one_ring_fs;
        vertex_one_ring_face_sides(G, v2fs.row(iV), one_ring_fs);
        double angle_sum = 0;
        for (const Vector2i& fs : one_ring_fs) {
            angle_sum += opposite_corner_angle(l, next(fs));
        }
        return angle_sum;
    };

    int nF = F.rows();
    int nV = v2fs.rows();

    // store local information
    Matrix<int, 3, 2> boundary_fs_twin;
    array<double, 3> radial_lengths, boundary_lengths, ang_coords, ang_scale;
    array<int, 3> faces{f, nF, nF+1}, vertices{F(f, 0), F(f, 1), F(f, 2)};
    for (int iS = 0; iS < 3; iS++) {
        boundary_fs_twin.row(iS) = twin(G, {f, iS});
        Vector3d edge_bary = b;
        edge_bary(iS) -= 1;
        radial_lengths[iS] = sqrt(barycentric_dot_product(edge_bary, edge_bary, l, f));
        boundary_lengths[iS] = l(f, iS);
        ang_coords[iS] = A(f, iS);

        int v = F(f, iS);
        ang_scale[iS] = (is_boundary_vertex(G, v2fs, v) ? M_PI : 2. * M_PI)
                        / vertex_angle_sum(F(f, iS));
    }

    // update the barycentric coordinates of all existing barycentric coordinates on face f
    // modify BC array in place
    array<vector<int>, 3> resident_points;
    for (int p : F2V[f]) {
        Vector3d x = BC.row(p);

        // barycentric coordinates in subface 0 are given by
        // {x[0] - x[2] / b[2] b[0], x[1] - x[2] / b[2] x[1], x[2] / x[1]}
        // pick subface where all are valid
        Vector3d c0{x[0] - x[2] / b[2] * b[0],
                    x[1] - x[2] / b[2] * b[1], x[2] / b[2]};
        Vector3d c1{x[1] - x[0] / b[0] * b[1],
                    x[2] - x[0] / b[0] * b[2], x[0] / b[0]};
        Vector3d c2{x[2] - x[1] / b[1] * b[2],
                    x[0] - x[1] / b[1] * b[0], x[1] / b[1]};

        if (!inside_triangle(c0) && !inside_triangle(c1) &&
            !inside_triangle(c2)) {
            cout << "subface barycentric coordinates for point " << x << ": " << c0.transpose()
                 << ", " << c1.transpose() << ", " << c2.transpose()<< endl;
            // TODO: round to closest triangle instead of throwing error
            throw runtime_error("point not in any subtriangle. maybe round?");
        }

        if (inside_triangle(c0)) {
            resident_points[0].push_back(p);
            BC.row(p) = c0;
        } else if (inside_triangle(c1)) {
            resident_points[1].push_back(p);
            BC.row(p) = c1;
        } else { // must be inside triangle c2
            resident_points[2].push_back(p);
            BC.row(p) = c2;
        }
    }

    // create two new faces and one new vertex
    // use conservativeResize to resize without losing old data
    F.conservativeResize(nF+2, 3);
    G.conservativeResize(nF+2, 6);
    l.conservativeResize(nF+2, 3);
    A.conservativeResize(nF+2, 3);
    v2fs.conservativeResize(nV+1, 2);
    int v = nV; // new vertex index is nV
    // new face indices are nF and nF+1

    // update face vertex adjacency
    F(f, 0) = vertices[0];
    F(f, 1) = vertices[1];
    F(f, 2) = v;

    F(nF, 0) = vertices[1];
    F(nF, 1) = vertices[2];
    F(nF, 2) = v;

    F(nF+1, 0) = vertices[2];
    F(nF+1, 1) = vertices[0];
    F(nF+1, 2) = v;

    // relabel boundary twin (so that self edges can be handled properly)
    Matrix<int, 3, 2> boundary_fs;
    boundary_fs.row(0) << f,    0;
    boundary_fs.row(1) << nF,   0;
    boundary_fs.row(2) << nF+1, 0;
    relabel(boundary_fs, f, boundary_fs_twin);

    // update gluing map
    for (int iS = 0; iS < 3; iS++)
      glue_face_sides(boundary_fs.row(iS), boundary_fs_twin.row(iS), G);
    glue_face_sides({f,    1}, {nF,   2}, G);
    glue_face_sides({nF,   1}, {nF+1, 2}, G);
    glue_face_sides({nF+1, 1}, {f,    2}, G);

    // update edge lengths
    l(f, 0) = boundary_lengths[0];
    l(f, 1) = radial_lengths[1];
    l(f, 2) = radial_lengths[0];

    l(nF, 0) = boundary_lengths[1];
    l(nF, 1) = radial_lengths[2];
    l(nF, 2) = radial_lengths[1];

    l(nF+1, 0) = boundary_lengths[2];
    l(nF+1, 1) = radial_lengths[0];
    l(nF+1, 2) = radial_lengths[2];

    // update angular coordinates
    //              v2
    //             / . \
    //            /  .  \
    //           /   .   \
    //      s2  /    .    \  s1
    //         /     .     \
    //        /      .      \
    //       / nF+1  b   nF  \
    //      /     .     .     \
    //     /   .     f     .   \
    //    / .                 . \
    //  v0  --------------------  v1
    //              s0
    A(f,    0) = ang_coords[0];
    A(nF,   0) = ang_coords[1];
    A(nF+1, 0) = ang_coords[2];

    A(f,    1) = ang_coords[1] + opposite_corner_angle(l, {nF,   1}) * ang_scale[1];
    A(nF,   1) = ang_coords[2] + opposite_corner_angle(l, {nF+1, 1}) * ang_scale[2];
    A(nF+1, 1) = ang_coords[0] + opposite_corner_angle(l, {f,    1}) * ang_scale[0];

    A(f,    2) = 0;
    A(nF,   2) = A(f,  2) + opposite_corner_angle(l, {f,  0});
    A(nF+1, 2) = A(nF, 2) + opposite_corner_angle(l, {nF, 0});

    // assign v2fs for new vertex, and ensure that v2fs always returns face side
    // with smallest angular coord for the other face vertices
    v2fs.row(v) << f, 2;
    v2fs.row(vertices[0]) = get_smallest_angular_coordinate(G, A, {f,    0});
    v2fs.row(vertices[1]) = get_smallest_angular_coordinate(G, A, {nF,   0});
    v2fs.row(vertices[2]) = get_smallest_angular_coordinate(G, A, {nF+1, 0});

    // update F2V lists
    F2V[f] = resident_points[0];
    F2V.push_back(resident_points[1]);
    F2V.push_back(resident_points[2]);

    // cout << " ----- inserted vertex " << v << " which is incident on faces "
    //      << f << ", " << nF << ", " << nF+1 << endl;

    return v;
}

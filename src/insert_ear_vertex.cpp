#include "insert_ear_vertex.h"

#include <array>

#include "next.h"
#include "twin.h"
#include "glue_face_sides.h"
#include "is_same_face_side.h"
#include "opposite_corner_angle.h"
#include "get_smallest_angular_coordinate.h"

int insert_ear_vertex(
                      const Eigen::Vector2i & fs,
                      double t,
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

    if (is_same_face_side(fs, GHOST_FACE_SIDE)) {
        throw std::invalid_argument("tried to insert a vertex onto a nonexistent"
                                    "boundary edge");
    }

    if (t < 0 || t > 1) {
        throw std::invalid_argument("vertex insertion point is outside of edge");
    }

    int nF = F.rows();
    int nV = v2fs.rows();

    // store local information
    array<int, 2> vertices{F(fs(0), fs(1)), F(fs(0), (fs(1)+1)%3)};
    double fs_len = l(fs(0), fs(1));
    array<double, 2> ang_coords{A(fs(0), fs(1)),
                                A(fs(0), (fs(1)+1)%3)
                                + opposite_corner_angle(l, next(next(fs)))};

    // create one new face and one new vertex
    // use conservativeResize to resize without losing old data
    F.conservativeResize(nF+1, 3);
    G.conservativeResize(nF+1, 6);
    l.conservativeResize(nF+1, 3);
    A.conservativeResize(nF+1, 3);
    v2fs.conservativeResize(nV+1, 2);
    int v = nV; // new vertex index is nV
    int f = nF; // new face index is nF

    // update face vertex adjacency
    F(f, 0) = vertices[1];
    F(f, 1) = vertices[0];
    F(f, 2) = v;

    // update gluing map
    glue_face_sides({f, 0}, {fs(0), fs(1)}, G);
    glue_face_sides({f, 1}, GHOST_FACE_SIDE, G);
    glue_face_sides({f, 2}, GHOST_FACE_SIDE, G);

    // update edge lengths
    l(f, 0) = fs_len;
    l(f, 1) = t * fs_len;
    l(f, 2) = (1-t) * fs_len;

    // update angular coordinates
    A(f, 0) = ang_coords[1];
    A(f, 1) = ang_coords[0];
    A(f, 2) = 0;

    // assign v2fs for new vertex, and ensure that v2fs always returns face side
    // with smallest angular coord for the other face vertices
    v2fs.row(v) << f, 2;
    v2fs.row(vertices[0]) << f, 1;

    // update F2V list. Since we're adding an ear triangle, it contains no points
    F2V.push_back(vector<int>{});

    return v;
}

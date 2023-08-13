#ifndef INSERT_ANGULAR_COORDINATE
#define INSERT_ANGULAR_COORDINATE

#include <Eigen/Core>
#include <Eigen/Dense>

#include "twin.h"
#include "next.h"
#include "opposite_corner_angle.h"
#include "cw.h"
#include "ccw.h"
#include "get_face_side_angular_coordinate.h"
#include "is_same_face_side.h"
#include "global_variables.h"
#include "vertex_one_ring_face_sides.h"
#include "pi.h"

// insert a new face side (without angular coordinate) to the angular coordinate of a vertex

// Inputs
//     G: Fx6 glue map
//     l: Fx3 edge lengths
//     fs: face side to be inserted
//     A: Fx3 angular coordinate matrix
// Output
//     A: upated A
void insert_angular_coordinate(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs,
    Eigen::MatrixXd & A);
#endif

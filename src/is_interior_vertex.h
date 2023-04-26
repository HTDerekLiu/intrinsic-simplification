#ifndef IS_INTERIOR_VERTEX
#define IS_INTERIOR_VERTEX

#include <Eigen/Core>
#include <Eigen/Dense>

#include "vertex_one_ring_face_sides.h"

// check whether the face side is in the interior
// Inputs
//     G: glue map
//     v2fs: vertex to face side map
//     v: a vertex
// Outputs
//     true/false 
bool is_interior_vertex(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v);
#endif
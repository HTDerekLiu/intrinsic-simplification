#ifndef IS_BOUNDARY_VERTEX
#define IS_BOUNDARY_VERTEX

#include <Eigen/Core>
#include <Eigen/Dense>

#include "vertex_one_ring_face_sides.h"

// check whether the vertex is on the boundary
// Inputs
//     G: glue map
//     v2fs: vertex to face side map
//     v: a vertex
// Outputs
//     true/false 
bool is_boundary_vertex(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v);
#endif
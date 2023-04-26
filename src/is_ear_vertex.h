#ifndef IS_EAR_VERTEX
#define IS_EAR_VERTEX

#include <Eigen/Core>
#include <Eigen/Dense>

#include "vertex_one_ring_face_sides.h"
#include "is_interior_vertex.h"

// check whether the vertex is an ear vertex
// Inputs
//     G: glue map
//     v2fs: vertex to face side map
//     v: a vertex
// Outputs
//     true/false 
bool is_ear_vertex(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v);
#endif
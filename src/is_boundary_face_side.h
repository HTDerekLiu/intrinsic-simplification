#ifndef IS_BOUNDARY_FACE_SIDE
#define IS_BOUNDARY_FACE_SIDE

#include <Eigen/Core>
#include <Eigen/Dense>

#include "global_variables.h"
#include "twin.h"
#include "is_same_face_side.h"

// check whether the face side is on the boundary
// Inputs
//     G: glue map
//     fs: a face side
// Outputs
//     true/false 
bool is_boundary_face_side(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs);
#endif
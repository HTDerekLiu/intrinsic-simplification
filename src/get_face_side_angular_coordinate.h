#ifndef GET_FACE_SIDE_ANGULAR_COORDINATE
#define GET_FACE_SIDE_ANGULAR_COORDINATE

#include <Eigen/Core>
#include <Eigen/Dense>
#include "global_variables.h"
#include "is_same_face_side.h"

// Get the angular coordinate of a given face side f, s
// Inputs
//     A: angular coordinates
//     fs: face-side
// Outputs
//     angle: angular coord
double get_face_side_angular_coordinate(
    const Eigen::MatrixXd & A,
    const Eigen::Vector2i & fs);
#endif
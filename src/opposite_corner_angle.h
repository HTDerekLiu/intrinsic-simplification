#ifndef OPPOSITE_CORNER_ANGLE
#define OPPOSITE_CORNER_ANGLE

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h> 
#include <algorithm>

#include "global_variables.h"

// Computes triangle corner angle opposite the face-side fs.
// Input
//     l: a |F|x3 array of face-side edge lengths
//     fs: a tuple of face side
// Outputs
//     angle:  The corner angle opposite to (f,s), in radians
double opposite_corner_angle(
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs);
#endif
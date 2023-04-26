#ifndef GET_ANGULAR_COORDINATE_RANGE
#define GET_ANGULAR_COORDINATE_RANGE

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <math.h> 

#include "is_boundary_vertex.h"
#include "vertex_one_ring_face_sides.h"
#include "opposite_corner_angle.h"
#include "twin.h"
#include "next.h"

// given a vertex v, this method outputs the maximum valid value in the angular coordinate. For interior vertex, this is always 2pi. For boundary vertex, this is [0,2pi]
// Inputs:
//     G: glue map
//     l: edge lengths
//     v2fs: vertex to fs map
//     v: vertex index
// Outputs:
//     angle_range
double get_angular_coordinate_range(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::MatrixXi & v2fs,
    const int & v);
#endif
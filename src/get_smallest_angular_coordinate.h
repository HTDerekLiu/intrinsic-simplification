#ifndef GET_SMALLEST_ANGULAR_COORDINATE
#define GET_SMALLEST_ANGULAR_COORDINATE

#include <Eigen/Core>
#include <Eigen/Dense>

#include "vertex_one_ring_face_sides.h"

#include <vector>

// loop over the face-side of a vertex and get the one with the smallest angular coordinate 

// Inputs
// G: |F|x6 array of gluing map. 
// A: |F|x3 array of signpost angles
// fs: a tuple of a face side

// Outputs:
// fs_min: a face side starting with vertex(F,fs) with smallest angular coordinates
Eigen::Vector2i get_smallest_angular_coordinate(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & A,
    const Eigen::Vector2i & fs);
#endif
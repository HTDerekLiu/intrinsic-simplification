#ifndef BUILD_ANGULAR_COORDINATES
#define BUILD_ANGULAR_COORDINATES

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <math.h>

#include <igl/cumsum.h>

#include "vertex_one_ring_face_sides.h"
#include "opposite_corner_angle.h"


// Builds the angular coordinates for each vertex in the mesh. Note that for the boundary vertices, our implementation assumes zero geodesic curvature. Thus it will result in large distortion if this is not the case.

// Inputs
//     F: |F|x3 array of face list
//     G: |F|x6 array of gluing map. 
//     l: |F|x3 array of face-side lengths

// Outputs:
//     A: |F|x3 array of angular coordinates. A(f,s) gives you the angular coordinate of this face-side from its starting vertex
//     v2fs: |V|x2 array of face-sides. v2fs.row(v) outputs one face-side starting from this vertex 

void build_angular_coordinates(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs);

#endif
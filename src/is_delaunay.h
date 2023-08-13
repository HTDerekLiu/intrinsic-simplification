#ifndef IS_DELAUNAY
#define IS_DELAUNAY

#include <Eigen/Core>
#include <Eigen/Dense>

#include <queue>
#include <vector>
#include <math.h>

#include "is_boundary_face_side.h"
#include "global_variables.h"
#include "opposite_corner_angle.h"
#include "twin.h"
#include "pi.h"

/*
    Test if the edge given by face-side fs satisfies the intrinsic Delaunay property.
    
    Inputs
    G: |F|x3x2 gluing map G,
    l: |F|x3 array of face-side edge lengths
    f: A integer of a face
    s: A integer of a side

    Outputs
    True/False on whether the edge is Delaunay
*/
bool is_delaunay(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs);
#endif

#ifndef UPDATE_ANGULAR_COORDINATE
#define UPDATE_ANGULAR_COORDINATE

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <math.h>

#include <igl/cumsum.h>

#include "is_boundary_vertex.h"
#include "vertex_one_ring_face_sides.h"
#include "opposite_corner_angle.h"
#include "twin.h"
#include "next.h"
#include "pi.h"

// this function updates the angular coordinate of a vertex with the same connecivity. This function will keep the angular coordinate to have the same origin as the original `A`
void update_angular_coordinate(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::MatrixXd & A);
#endif

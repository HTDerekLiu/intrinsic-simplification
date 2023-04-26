#ifndef FLATTEN_INTERIOR_VERTEX_AND_COST
#define FLATTEN_INTERIOR_VERTEX_AND_COST

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h> 
#include <set>
#include <tuple>
#include <vector>
#include <unordered_set>

#include "is_interior_vertex.h"
#include "vertex_one_ring_face_sides.h"
#include "vertex_one_ring_unique_face_sides.h"
#include "opposite_corner_angle.h"
#include "vertex_one_ring_vertices.h"
#include "vertex_one_ring_faces.h"
#include "get_face_side_vertices.h"
#include "gaussian_curvature_at_vertex.h"
#include "polar_to_cartesian.h"
#include "cetm_flatten_interior_vertex.h"
#include "global_variables.h"
#include "update_angular_coordinate.h"
#include "is_one_ring_fs_valid.h"

// pre flattening a vertex and compute the flattening cost. Note that this will not change the input arguments, it only computes the decimation cost
double pre_flatten_interior_vertex_and_cost(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    Eigen::MatrixXd & l, // we will touch this, but we won't change its values
    const Eigen::MatrixXd & A,
    const Eigen::MatrixXi & v2fs,
    const Eigen::MatrixXd & T,
    const Eigen::VectorXd & K,
    const int & v_flatten);

double flatten_interior_vertex_and_cost(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const std::vector<std::vector<int>> & F2V,
    const int & v_flatten,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXd & BC,
    Eigen::MatrixXd & T,
    Eigen::VectorXd & K);

#endif
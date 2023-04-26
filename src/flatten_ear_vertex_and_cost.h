#ifndef FLATTEN_EAR_VERTEX_AND_COST
#define FLATTEN_EAR_VERTEX_AND_COST

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h> 
#include <set>
#include <tuple>
#include <vector>

#include "flatten_boundary_vertex_and_cost.h"
#include "is_ear_vertex.h"
#include "vertex_one_ring_face_sides.h"
#include "flip_edge.h"
#include "get_face_side_vertices.h"
#include "glue_face_sides.h"
#include "flip_edge.h"
#include "flip_edge_cw.h"
#include "global_variables.h"
#include "is_same_face_side.h"

// pre flattening a vertex and compute the flattening cost. Note that this will not change the input arguments, it only computes the decimation cost
double pre_flatten_ear_vertex_and_cost(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l, 
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    const Eigen::MatrixXd & T,
    const Eigen::VectorXd & K,
    const int & v_flatten);

double flatten_ear_vertex_and_cost(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXi & v2fs,
    std::vector<std::vector<int>> & F2V,
    const int & v_flatten,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXd & BC,
    Eigen::MatrixXd & T,
    Eigen::VectorXd & K);

#endif
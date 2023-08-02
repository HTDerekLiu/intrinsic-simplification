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

/*
pre flattening an ear vertex and compute the flattening cost. Note that this will not change the input arguments, it only computes the decimation cost

Inputs
    F: |F|x3 vertex-face adjacency list
    G: |F|x6 glue map
    l: |F|x3 edge lengths for each face side
    A: |F|x3 angular coordinate for each face side
    v2fs: |v|x2 where v2fs.row(i) returns a face side for vertex i
    T: |V|x(num_quantitiesx3) transport cost information, where T[v,0] is the mass, T[v,2:3] is the vector to karcher min
    K: Gaussian curvature at each vertex (cache some computation for efficiency reasons)
    v_flatten: vertex index that needs to be flattened

Output
    cost of flattening v_flatten
*/
double pre_flatten_ear_vertex_and_cost(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l, 
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    const Eigen::MatrixXd & T,
    const Eigen::VectorXd & K,
    const int & v_flatten);

/*
flattening an ear vertex and compute the flattening cost. The difference with the above function is that all the intrinsic information will get updated in place

Inputs
    F: |F|x3 vertex-face adjacency list
    G: |F|x6 glue map
    l: |F|x3 edge lengths for each face side
    A: |F|x3 angular coordinate for each face side
    v2fs: |v|x2 where v2fs.row(i) returns a face side for vertex i
    T: |V|x(num_quantitiesx3) transport cost information, where T[v,0] is the mass, T[v,2:3] is the vector to karcher min
    K: Gaussian curvature at each vertex (cache some computation for efficiency reasons)
    v_flatten: vertex index that needs to be flattened

Output
    cost of flattening v_flatten
    l, A, BC, T, K are changed in place
*/
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
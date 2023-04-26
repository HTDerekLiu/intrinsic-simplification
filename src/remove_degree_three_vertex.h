#ifndef REMOVE_DEGREE_THREE_VERTEX
#define REMOVE_DEGREE_THREE_VERTEX

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h> 
#include <set>
#include <stdexcept>

#include "is_interior_vertex.h"
#include "vertex_one_ring_face_sides.h"
#include "get_face_side_vertices.h"
#include "twin.h"
#include "next.h"
#include "roll1d.h"
#include "global_variables.h"
#include "glue_face_sides.h"
#include "get_smallest_angular_coordinate.h"
#include "flip_to_delaunay_local.h"
#include "is_one_ring_fs_valid.h"

/*
    This function removes a degree three vertex

    Inputs
    v: the vertex index we try to remove

    V: |V|x3 array of vertex locations
    F: |F|x3 array of face list
    G: |F|x3x2 array of gluing map 
    l: |F|x3 array of edge-lengths array
    A: |F|x3 array of signpost angles
    v2fs: |V|x2 array of face-sides
    BC: |BC|x3 array of barycentric coordinates 
    F2V: |F| length list of lists, where F2V[f] gives you indices in BC.
    T: |V|xch*3 array of tangent vectors

    Outputs
    All the variables are updated in-place
*/
void remove_degree_three_vertex(
    const int & v,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V,
    Eigen::MatrixXd & T);

#endif
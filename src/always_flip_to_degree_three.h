#ifndef ALWAYS_FLIP_TO_DEGREE_THREE
#define ALWAYS_FLIP_TO_DEGREE_THREE

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h> 
#include <set>
#include <stdexcept>

#include "is_interior_vertex.h"
#include "vertex_one_ring_face_sides.h"
#include "opposite_corner_angle.h"
#include "vertex_one_ring_vertices.h"
#include "vertex_one_ring_faces.h"
#include "get_face_side_vertices.h"
#include "twin.h"
#include "next.h"
#include "is_diamond_convex.h"
#include "flip_edge.h"
#include "remove_vector_element.h"
#include "global_variables.h"
#include "gaussian_curvature_at_vertex.h"
#include "is_one_ring_fs_valid.h"

/*
    Given a vertex v ***which can be flipped to degree three***, this function flips edges so that it has degree three. If this vertex has zero Gaussian curvature, then it can always be flipped to degree-3

    Inputs:
    F: |F|x3 vertex-face adjacency list F
    G: |F|x3x2 gluing map G
    El: |F|x3 edge-lengths array, giving the length of each face-side
    A: |F|x3 array of signpost angles
    v2fs: |V|x2 array of face-sides. V2FS[v] outputs one face-side starting from this vertex 
    BC: |BC|x3 array of barycentric coordinates whose corresponding faces are stored in F2V implicitly
    F2V: |F| length list of lists, where F2V[f] gives you a list of indices in BC. For example, if F2V[f] = [v], then BC[v,:] corresponds to the barycentric coordinates in F[f,:]
    v: a vertex index

    Outputs 
    all the intrinsic information changed in place
*/
void always_flip_to_degree_three(
    const int & v,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V);

#endif
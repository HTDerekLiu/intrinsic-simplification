#ifndef ALWAYS_FLIP_TO_EAR
#define ALWAYS_FLIP_TO_EAR

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h> 
#include <set>
#include <stdexcept>
#include <iterator>
#include <vector>

#include "is_boundary_vertex.h"
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
    Given a boundary vertex v which CAN BE FLIPPED, this function flips edges so that the vertex becomes an ear

    Inputs:
    F: |F|x3 vertex-face adjacency list F
    G: |F|x3x2 gluing map G
    El: |F|x3 edge-lengths array, giving the length of each face-side
    A: |F|x3 array of signpost angles
    v2fs: |V|x2 array of face-sides. V2FS[v] outputs one face-side starting from this vertex 
    v: a vertex index

    Additional inputs:
    BC: |BC|x3 array of barycentric coordinates whose corresponding faces are stored in F2V implicitly
    F2V: |F| length list of lists, where F2V[f] gives you a list of indices in BC. For example, if F2V[f] = [v], then BC[v,:] corresponds to the barycentric coordinates in F[f,:]
    ophist: operation history list, appended-to in-place

    Outputs 
    Ture/False indicates whether the edge flip successfully flip to an ear
*/
void always_flip_to_ear(
    const int & v,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V);

#endif
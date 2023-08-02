#ifndef FLIP_TO_DELAUNAY
#define FLIP_TO_DELAUNAY

#include <Eigen/Core>
#include <Eigen/Dense>

#include <queue>
#include <vector>
#include <iostream>

#include "is_diamond_convex.h"
#include "flip_edge.h"
#include "is_delaunay.h"
#include "is_boundary_face_side.h"
#include "twin.h"
#include "next.h"

/*
Flip edges in the triangulation until it satisifes the intrinsic Delaunay criterion.

Inputs
    F: |F|x3 array of face list
    G: |F|x3x2 array of gluing map 
    l: |F|x3 array of edge-lengths array
    A: |F|x3 array of signpost angles
    V2FS: |V|x2 array of face-sides. V2FS[v] outputs one face-side starting from this vertex 

Optional inputs:
    BC: |BC|x3 array of barycentric coordinates whose corresponding faces are stored in F2V implicitly
    F2V: |F| length list of lists, where F2V[f] gives you a list of indices in BC. For example, if F2V[f] = [v], then BC[v,:] corresponds to the barycentric coordinates in F[f,:]
    ophist: operation history list, appended-to in-place

Outputs
    changed in place
*/
void flip_to_delaunay(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V);

void flip_to_delaunay(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs);
#endif
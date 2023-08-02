#ifndef FLIP_EDGE_CW
#define FLIP_EDGE_CW

#include <Eigen/Core>
#include <Eigen/Dense>

#include <assert.h>
#include <map>
#include <stdexcept>
#include <vector>

#include "is_boundary_face_side.h"
#include "is_diamond_convex.h"
#include "compute_barycentric_robust.h"
#include "flatten_diamond_mesh.h"
#include "get_smallest_angular_coordinate.h"
#include "insert_angular_coordinate.h"
#include "vertex_one_ring_face_sides.h"
#include "glue_face_sides.h"
#include "roll1d.h"
#include "diagonal_length.h"
#include "get_face_side_angular_coordinate.h"
#include "twin.h"
#include "next.h"
#include "is_same_face_side.h"

/*
Performs an intrinsic edge flip on the edge given by face-side (f, s). This edge flip rotate the (f,s) clockwise.

Inputs
    s0: A face side to be flipped
    F: |F|x3 vertex-face adjacency list F
    G: |F|x3x2 gluing map G
    El: |F|x3 edge-lengths array
    A: |F|x3 array of angular coordinates of each face side
    v2fs: |V|x2 array of face-sides. v2fs[v] outputs one face-side starting from this vertex 

Optional inputs:
    BC: |BC|x3 array of barycentric coordinates
    F2V: |F| length list of lists, where F2V[f] gives you a list of indices in BC.

Outputs 
    everything updated in place besides s0
    fs_dict: store the map between old fs (as std::pair) and new fs (as std::pair). Note that this is different from the usual represnetation Vector2i

TODO: add output 
ophist: a object of operation info for undoing the flip

Notation
         before          after  
           v2             v2
          /  \           / | \  
         /    \         /  |  \    
        /s2  s1\       /s2 | s1\ 
       /        \     /    |    \ 
      /    s0    \   /     |     \  
     v0 -------- v1 v0  t0 | s0  v1
      \    t0    /   \     |     /
       \        /     \    |    /  
        \t1  t2/       \t1 | t2/      
         \    /         \  |  /    
          \  /           \ | /
           v3             v3
*/
void flip_edge_cw(
    const Eigen::Vector2i & s0,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V,
    std::map<std::pair<int, int>, std::pair<int, int>> & fs_dict);

void flip_edge_cw(
    const Eigen::Vector2i & s0,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V);

void flip_edge_cw(
    const Eigen::Vector2i & s0,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs);
#endif
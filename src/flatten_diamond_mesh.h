#ifndef FLATTEN_DIAMOND_MESH
#define FLATTEN_DIAMOND_MESH

#include <Eigen/Core>
#include <Eigen/Dense>

#include "next.h"
#include "twin.h"
#include "opposite_corner_angle.h"
#include "roll1d.h"

/*
flatten the diamond of a face side into 2D 

Inputs
    G: glue map
    l: edge length matrix
    fs: face side

Optional inputs
    is_ccw: true/false to determine whether edge flip is ccw or cw

Outputs
    U: 4x2 UV coordinates
    F: face indices before flip

Optional outputs:
    F_flip: face indices after flip (the "FF" in the figure)

Notation of ccw flip
         before          after  
           u2              u2
          /  \            / | \  
         /    \          /  |  \    
        / F[0] \        / fs|   \ 
       /        \      /    |    \ 
      /    fs    \    /     |     \  
     u0 -------- u1 u0 FF[0]|FF[1] u1
      \          /    \     |     /
       \        /      \    |    /  
        \ F[1] /        \   |   /      
         \    /          \  |  /    
          \  /            \ | /
           u3              u3
*/

void flatten_diamond_mesh(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs,
    Eigen::MatrixXd & U,
    Eigen::MatrixXi & F);

void flatten_diamond_mesh(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs,
    const bool & is_ccw, // is "ccw" or "cw" flip direction
    Eigen::MatrixXd & U,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & F_flip);
#endif
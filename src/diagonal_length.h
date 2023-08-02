#ifndef DIAGONAL_LENGTH
#define DIAGONAL_LENGTH

#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h> 
#include <iostream> 

#include "next.h"
#include "twin.h"
#include "opposite_corner_angle.h"

/*
Computes the length of the opposite diagonal of the diamond formed by the triangle containing fs, and the neighboring triangle adjacent to fs.

Inputs
    G: |F|x3x2 gluing map
    l: |F|x3 array of face-side edge lengths
    fs: A face-side (f,s)

Outputs:
    The diagonal length

Note:
         / \
        / | \  
       /  |  \
      p   d   u
     /    |    \
    /     fs  a \
    ------------- 
    \     |   b /
     \    |    /
      q   |   v 
       \  |  /
        \ | /
         \ /
*/
double diagonal_length(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs);
#endif
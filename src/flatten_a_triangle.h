#ifndef FLATTEN_A_TRIANGLE
#define FLATTEN_A_TRIANGLE

#include <Eigen/Core>
#include <Eigen/Dense>

/*
given three edge lengths of a triangle, this function computes "a" 2D embedding of that triangle

Inputs:
    lij, ljk, lkl: lengths of the triangle edges

Output:
    UV: 3x2 UV locations of the triangle 

Notes:
            v0 
           / \
          /   \  
         /     \
       l[0]     l[2]
       /         \
      /           \
     v1 ---l[1]--- v2
    
    v0 = [0, 0]
    v1 = [l0, 0] 
*/
void flatten_a_triangle(
    const double & lij,
    const double & ljk,
    const double & lkl,
    Eigen::MatrixXd & UV);
#endif
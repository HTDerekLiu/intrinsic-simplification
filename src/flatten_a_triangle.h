#ifndef FLATTEN_A_TRIANGLE
#define FLATTEN_A_TRIANGLE

#include <Eigen/Core>
#include <Eigen/Dense>

/*
    given three edge lengths of a triangle, this function computes "a" 2D embedding of that triangle

    Inputs:
    l: (3,) array of the length of each edge of a triangle

    Output:
    v0: (2,) array fo the location of a tri corner point of the triangle (we set v0 == [0,0])
    v1: (2,) array of the location of the next tri corner point, such that len(v1,v0) = l0 and we let v1 = [l0, 0]
    v2: (2,) array of the location of the last tri corner point, which is uniquely determined by the edge length

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
#ifndef COTAN_LAPLACIAN
#define COTAN_LAPLACIAN

#include <Eigen/Sparse>

#include <vector>
#include <math.h>
#include <iostream>

#include "opposite_corner_angle.h"

/*
Compute cotangent weighed Laplcain

Inputs: 
  F: |F|x3 vertex-face adjacency list
  l: |F|x3 edge lengths for each face side
Outputs:
  L: |V|x|V| cotan laplace matrix
*/
void cotan_Laplacian(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & l,
    Eigen::SparseMatrix<double> & L);

#endif
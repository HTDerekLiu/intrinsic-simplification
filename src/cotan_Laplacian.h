#ifndef COTAN_LAPLACIAN
#define COTAN_LAPLACIAN

#include <Eigen/Sparse>

#include <vector>
#include <math.h>
#include <iostream>

#include "opposite_corner_angle.h"

// Compute cotangent weighed Laplcain
//
// Inputs: 
//   l face side edge lengths
//   F face indices
// Outputs:
//   L cotan laplace matrix
void cotan_Laplacian(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & l,
    Eigen::SparseMatrix<double> & L);

#endif
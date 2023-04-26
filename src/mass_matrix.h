#ifndef MASS_MATRIX
#define MASS_MATRIX

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "vertex_areas.h"

void mass_matrix(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & l,
    Eigen::SparseMatrix<double> & M);

#endif
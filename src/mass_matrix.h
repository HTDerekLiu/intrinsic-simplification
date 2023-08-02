#ifndef MASS_MATRIX
#define MASS_MATRIX

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "vertex_areas.h"

/*
Computes the lumped mass matrix of vertex areas

Inputs:
    F: |F|x3 face list
    l: |F|x3 face side lengths
Outputs:
    M: |V|x|V| sparse diagonal matrix of barycentric vertex areas
*/
void mass_matrix(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & l,
    Eigen::SparseMatrix<double> & M);

#endif
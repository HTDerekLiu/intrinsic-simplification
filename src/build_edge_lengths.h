#ifndef BUILD_EDGE_LENGTHS
#define BUILD_EDGE_LENGTHS

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>

/*
Build edge lengths for each face-side
Inputs
    V vertex locations
    F face indices 
Outputs
    l: |F|x3 matrix of face-side (half edge) lengths
*/
void build_edge_lengths(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & l);
#endif
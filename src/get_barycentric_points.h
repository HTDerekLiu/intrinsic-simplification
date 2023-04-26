#ifndef GET_BARYCENTRIC_POINTS
#define GET_BARYCENTRIC_POINTS

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

// turn barycentric points (BF, F2V) into points in R3

// Inputs:
//     V: |V|x3 array of vertex locations
//     F: |F|x3 array of face list
//     BC: nx3 barycentric coordinates
//     F2V: F2V[f] gives us a list of indices to BC on face f

// Outputs:
//     P: nx3 barycentric points in R3
void get_barycentric_points(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & BC,
    const std::vector<std::vector<int>> & F2V,
    Eigen::MatrixXd & P);

void get_barycentric_points_redundant(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & BC,
    const std::vector<std::vector<int>> & F2V,
    Eigen::MatrixXd & P);
#endif
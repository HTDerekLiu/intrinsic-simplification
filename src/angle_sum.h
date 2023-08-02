#ifndef ANGLE_SUM
#define ANGLE_SUM

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

#include <opposite_corner_angle.h>
#include <next.h>
#include <global_variables.h>

/*
ANGLE_SUM computes the angle sum for all vertices

Inputs:
    F: |F|x3 vertex-face adjacency list F
    l: |F|x3 edge lengths for each face side

Outputs 
    ang_sum: |V| vector of all vertex angle sums
*/
void angle_sum(
    const Eigen::MatrixXi &F, 
    const Eigen::MatrixXd &l, 
    Eigen::VectorXd &ang_sum);

/*
this is the helper function for the above angle_sum
*/
void angle_sum(
    const int & nV,
    const Eigen::MatrixXi &F, 
    const Eigen::MatrixXd &l, 
    Eigen::VectorXd &ang_sum);

#endif

#ifndef ANGLE_SUM
#define ANGLE_SUM

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

#include <opposite_corner_angle.h>
#include <next.h>
#include <global_variables.h>

void angle_sum(
    const Eigen::MatrixXi &F, 
    const Eigen::MatrixXd &l, 
    Eigen::VectorXd &ang_sum);

void angle_sum(
    const int & nV,
    const Eigen::MatrixXi &F, 
    const Eigen::MatrixXd &l, 
    Eigen::VectorXd &ang_sum);

#endif

#ifndef ROTATION_MATRIX_AROUND_AXIS
#define ROTATION_MATRIX_AROUND_AXIS

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h> 

/*
generate a 3x3 rotation matrix to rotate a vector in R3 aound an axis for angle theta

Inputs
theta: a scalar angle
u: (3,) numpy array of rotation axis

Outputs
R: 3x3 rotation matrix 
*/

void rotation_matrix_around_axis(
    const double & theta,
    const Eigen::VectorXd & u,
    Eigen::MatrixXd & R);

void rotation_matrix_around_axis(
    const double & theta,
    const Eigen::Vector3d & u,
    Eigen::Matrix3d & R);

#endif
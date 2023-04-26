#include "rotation_matrix_around_axis.h"

void rotation_matrix_around_axis(
    const double & theta,
    const Eigen::VectorXd & u,
    Eigen::MatrixXd & R)
{
    double inv_u_norm = 1.0 / u.norm();
    double ux = u(0) * inv_u_norm;
    double uy = u(1) * inv_u_norm;
    double uz = u(2) * inv_u_norm;

    double c = cos(theta);
    double s = sin(theta);

    R.resize(3,3);
    R << c + ux*ux*(1-c)   , ux*uy*(1-c) - uz*s, ux*uz*(1-c) + uy*s,
         uy*ux*(1-c) + uz*s, c+uy*uy*(1-c)     , uy*uz*(1-c)-ux*s  ,
         uz*ux*(1-c) - uy*s, uz*uy*(1-c) + ux*s, c + uz*uz*(1-c);
}

void rotation_matrix_around_axis(
    const double & theta,
    const Eigen::Vector3d & u,
    Eigen::Matrix3d & R)
{
    double inv_u_norm = 1.0 / u.norm();
    double ux = u(0) * inv_u_norm;
    double uy = u(1) * inv_u_norm;
    double uz = u(2) * inv_u_norm;

    double c = cos(theta);
    double s = sin(theta);

    // R.resize(3,3);
    R << c + ux*ux*(1-c)   , ux*uy*(1-c) - uz*s, ux*uz*(1-c) + uy*s,
         uy*ux*(1-c) + uz*s, c+uy*uy*(1-c)     , uy*uz*(1-c)-ux*s  ,
         uz*ux*(1-c) - uy*s, uz*uy*(1-c) + ux*s, c + uz*uz*(1-c);
}
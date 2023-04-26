#include "flatten_a_triangle.h"
void flatten_a_triangle(
    const double & lij,
    const double & ljk,
    const double & lki,
    Eigen::MatrixXd & UV)
{
    UV.resize(3,2);
    UV.row(0) << 0, 0;
    UV.row(1) << lij, 0;
    double ang_kij = acos((lij*lij + lki*lki - ljk*ljk) / (2*lij*lki));
    UV.row(2) << lki * cos(ang_kij), lki * sin(ang_kij);
}
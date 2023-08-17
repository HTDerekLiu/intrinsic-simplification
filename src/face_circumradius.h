#ifndef FACE_CIRCUMRADIUS
#define FACE_CIRCUMRADIUS

#include <Eigen/Core>
#include <Eigen/Dense>

/*
compute the circumradius of a face

Inputs:
    l: |F|x3 edge lengths for each face side
    iF: index of face to compute circumradius for

Outputs:
    returns circumradius of face iF
*/

double face_circumradius(
    const Eigen::MatrixXd & l,
    int iF);

#endif
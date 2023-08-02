#ifndef FACE_AREAS
#define FACE_AREAS

#include <Eigen/Core>
#include <Eigen/Dense>

/*
compute area for each face

Inputs:
    l: |F|x3 edge lengths for each face side

Outputs:
    FA: |F| vector of face areas

*/

void face_areas(
    const Eigen::MatrixXd & l,
    Eigen::VectorXd & FA);

#endif
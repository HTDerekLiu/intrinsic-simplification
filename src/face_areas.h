#ifndef FACE_AREAS
#define FACE_AREAS

#include <Eigen/Core>
#include <Eigen/Dense>

// compute face area for each face

void face_areas(
    const Eigen::MatrixXd & l,
    Eigen::VectorXd & FA);

#endif
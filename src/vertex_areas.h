#ifndef VERTEX_AREAS
#define VERTEX_AREAS

#include <Eigen/Core>
#include <Eigen/Dense>

#include "face_areas.h"

// compute face area for each face

void vertex_areas(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & l,
    Eigen::VectorXd & VA);

#endif
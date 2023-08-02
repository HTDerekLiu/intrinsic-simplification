#ifndef VERTEX_AREAS
#define VERTEX_AREAS

#include <Eigen/Core>
#include <Eigen/Dense>

#include "face_areas.h"

// compute vertex area for each vertex
// Inputs
// F: |F|x3 face list
// l: |F|x3 face side lengths

// Outputs
// VA: |V| vector of vertex areas

void vertex_areas(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & l,
    Eigen::VectorXd & VA);

#endif
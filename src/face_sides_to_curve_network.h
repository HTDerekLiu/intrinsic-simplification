#ifndef FACE_SIDES_TO_CURVE_NETWORK
#define FACE_SIDES_TO_CURVE_NETWORK

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

#include "get_face_side_vertices.h"

/*
given a list of face sides, this funvtio creates a curve network for visualization.
*/
void face_sides_to_curve_network(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const std::vector<Eigen::Vector2i> & fs_list,
    Eigen::MatrixXd & nodes,
    Eigen::MatrixXi & edges);

#endif
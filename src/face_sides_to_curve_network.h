#ifndef FACE_SIDES_TO_CURVE_NETWORK
#define FACE_SIDES_TO_CURVE_NETWORK

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

#include "get_face_side_vertices.h"

/*
given a list of face sides, this funvtio creates a curve network for visualization.

Inputs:
    V: |V|x3 vertex list
    F: |F|x3 vertex-face adjacency list
    fs_list: vector of Vector2i, it contains a list of face sides

Inputs
    nodes: |P|x3 node location list
    edges: |P|x2 node edge list
*/
void face_sides_to_curve_network(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const std::vector<Eigen::Vector2i> & fs_list,
    Eigen::MatrixXd & nodes,
    Eigen::MatrixXi & edges);

#endif
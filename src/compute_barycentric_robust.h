#ifndef COMPUTE_BARYCENTRIC_ROBUST
#define COMPUTE_BARYCENTRIC_ROBUST

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <algorithm>
#include "global_variables.h"

/*
Givne a 2d point "p" and the vertex locations of a triangle (a,b,c), this function outputs the barycentric coordinates

Inputs:
    p: 2D location of query point
    a,b,c: 2D corner locations of a triangle

Outputs:
    b: barycentric coordinates

Optional outputs:
    is_degenerated: whether triangle is degenerated
    is_inside: whether p is inside the triangle
    distance_to_valid: distance of the raw b to a valid barycentric coordinate

Reference:
"Barycentric coordinates computation in homogeneous coordinates" by Vaclav Skala
*/

Eigen::Vector3d compute_barycentric_robust(
    const Eigen::Vector2d & p,
    const Eigen::Vector2d & a,
    const Eigen::Vector2d & b,
    const Eigen::Vector2d & c);

Eigen::Vector3d compute_barycentric_robust(
    const Eigen::Vector2d & p,
    const Eigen::Vector2d & a,
    const Eigen::Vector2d & b,
    const Eigen::Vector2d & c,
    const double & degenerate_threshold);

Eigen::Vector3d compute_barycentric_robust(
    const Eigen::Vector2d & p,
    const Eigen::Vector2d & a,
    const Eigen::Vector2d & b,
    const Eigen::Vector2d & c,
    bool & is_degenerated,
    bool & is_inside,
    double & distance_to_valid);

Eigen::Vector3d compute_barycentric_robust(
    const Eigen::Vector2d & p,
    const Eigen::Vector2d & a,
    const Eigen::Vector2d & b,
    const Eigen::Vector2d & c,
    const double & degenerate_threshold,
    bool & is_degenerated,
    bool & is_inside,
    double & distance_to_valid);
#endif
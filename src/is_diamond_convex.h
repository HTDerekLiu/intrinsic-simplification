#ifndef IS_DIAMOND_CONVEX
#define IS_DIAMOND_CONVEX

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

#include "global_variables.h"
#include "twin.h"
#include "next.h"
#include "opposite_corner_angle.h"
#include "pi.h"

// this function check the diamond of the face-side is convex or not. Specifically, check whether angle_v0 or angle_v1 is larger than np.pi
// Input
//     G: |F|x3x2 gluing map G
//     l: |F|x3 edge-lengths array, giving the length of each face-side
//     fs: A face side
// Output
//     True/False of whether this diamond is convex
// Notation
// fs is the face side from i>j, and the opposite vertex of fs is k. The opposite vertex of twin(fs) is l
bool is_diamond_convex(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::Vector2i & fs);
#endif

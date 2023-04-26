#ifndef IS_SAME_FACE_SIDE
#define IS_SAME_FACE_SIDE

#include <Eigen/Core>
#include <Eigen/Dense>

// Given two face-sides (fs0, fs1), returns true/false on whether fs1 is the same as fs0
// Inputs
//     fs0: a face side
//     fs1: another face side
// Outputs
//     return true/false
bool is_same_face_side(
    const Eigen::Vector2i & fs0,
    const Eigen::Vector2i & fs1);
#endif
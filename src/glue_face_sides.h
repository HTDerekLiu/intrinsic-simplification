#ifndef GLUE_FACE_SIDES
#define GLUE_FACE_SIDES

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include "global_variables.h"
#include "is_same_face_side.h"

// glue two face sides (fs0, fs1) in G
// Inputs
//     fs0: face side #0
//     fs1: face side #1
// Outputs
//     G: updated gluing map where fs0, fs1 are twin face sides
//
// Note:
// G.block(f,s*2,1,2) is the way to access the twin of (f,s)

void glue_face_sides(
    const Eigen::Vector2i & fs0,
    const Eigen::Vector2i & fs1,
    Eigen::MatrixXi & G);
#endif
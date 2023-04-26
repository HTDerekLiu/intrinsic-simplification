#ifndef IS_DIAMOND_LOOPED
#define IS_DIAMOND_LOOPED

#include <Eigen/Core>
#include <Eigen/Dense>

#include <set>

#include "twin.h"
#include "next.h"

// check whether the diamond is a looped (v0~v3 are not four distinct vertices)
// Inputs
//     F: face list
//     G: glue map
//     s0: a face side
// Outputs
//     true/false 
bool is_diamond_looped(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & s0);
#endif
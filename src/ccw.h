#ifndef CCW
#define CCW

#include <Eigen/Core>
#include <Eigen/Dense>

#include "next.h"
#include "twin.h"
#include "fast_mod.h"

// For a given face-side fs, returns the counter clock wise rotated face side
// Inputs
//     G: |F|x6 gluing map G,
//     fs: a face-side tuple (f,s)
// Outputs
//     fs_ccw = fs.next.next.twin
Eigen::Vector2i ccw(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs);

#endif
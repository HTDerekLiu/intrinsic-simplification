#ifndef CW
#define CW

#include <Eigen/Core>
#include <Eigen/Dense>

#include "next.h"
#include "twin.h"
// #include "fast_mod.h"
#include "global_variables.h"

/*
For a given face-side fs, returns the clock wise rotated face side
Inputs
    G: |F|x6 gluing map G,
    fs: a face-side tuple (f,s)
Outputs
    fs_cw = fs.twin.next
*/
Eigen::Vector2i cw(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs);

#endif
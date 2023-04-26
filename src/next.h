#ifndef NEXT
#define NEXT

#include <Eigen/Core>
#include <Eigen/Dense>
#include "fast_mod.h"
#include "global_variables.h"

// For a given face-side fs, returns the next face-side. 
// Inputs
//     fs: A face side (f,s)
// Outputs
//     The next face-side in the same triangle 
Eigen::Vector2i next(
    const Eigen::Vector2i & fs);

#endif
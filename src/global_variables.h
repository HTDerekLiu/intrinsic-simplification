#ifndef GLOBAL_VARIABLES
#define GLOBAL_VARIABLES

#include <Eigen/Core>
#include <Eigen/Dense>
#include <limits.h>
#include <vector>

namespace global_variables
{
    // define some useful constants used across the library
    const int GHOST_INDEX = INT_MIN;
    const double EPS = 1e-5;
    const double DOUBLE_INF = std::numeric_limits<double>::infinity();
    const double DOUBLE_NAN = std::numeric_limits<double>::quiet_NaN();
    const double SELF_EDGE_SCORE = 6.3; // edge flip score for self edges

    const std::vector<int> GHOST_FACE_SIDE_data = {INT_MIN,INT_MIN}; 
    const Eigen::Vector2i GHOST_FACE_SIDE(GHOST_FACE_SIDE_data.data());

}

#endif
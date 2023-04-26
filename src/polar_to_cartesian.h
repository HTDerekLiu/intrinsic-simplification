#ifndef POLAR_TO_CARTESIAN
#define POLAR_TO_CARTESIAN

#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h> 

// transform polar coordiante to cartiesian coordinate
void polar_to_cartesian(
    const double & rho,
    const double & phi,
    Eigen::Vector2d & xy);
#endif
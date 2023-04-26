#include "polar_to_cartesian.h"

void polar_to_cartesian(
    const double & rho,
    const double & phi,
    Eigen::Vector2d & xy)
{
    xy << rho*cos(phi), rho*sin(phi);
}
#include "barycentric_vector_to_cartesian.h"

std::complex<double> barycentric_vector_to_cartesian(
                                                     const Eigen::Vector3d & u,
                                                     const Eigen::MatrixXd & l,
                                                     const int & f
                                                     )
{

  double u_length = sqrt(abs(barycentric_dot_product(u, u, l, f)));
  Eigen::Vector3d x_axis{-1,1,0};

  double cosTheta = barycentric_dot_product(u, x_axis, l, f) / (u_length * l(f, 0));
  cosTheta = fmin(fmax(cosTheta, -1),1); // clamp to [-1, 1]

  // use sign of third barycentric coordinate to determine sign of angle
  double theta = copysign(acos(cosTheta), u(2));

  return std::polar(u_length, theta);
}

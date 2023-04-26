#include "barycentric_dot_product.h"

double barycentric_dot_product(
                               const Eigen::Vector3d & u,
                               const Eigen::Vector3d & v,
                               const Eigen::MatrixXd & l,
                               const int & f
                               )
{
  double l0 = l(f, 1); // length opposite vertex 0
  double l1 = l(f, 2); // length opposite vertex 1
  double l2 = l(f, 0); // length opposite vertex 2

  // from geometrycentral/surface/barycentric_vector.ipp
  double term1 = (u(1) * v(0) + u(0) * v(1)) * l2 * l2;
  double term2 = (u(2) * v(1) + u(1) * v(2)) * l0 * l0;
  double term3 = (u(2) * v(0) + u(0) * v(2)) * l1 * l1;
  return -0.5 * (term1 + term2 + term3);
}

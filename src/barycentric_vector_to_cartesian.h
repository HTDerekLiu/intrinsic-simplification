#ifndef BARYCENTRIC_VECTOR_TO_CARTESIAN
#define BARYCENTRIC_VECTOR_TO_CARTESIAN

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <set>
#include <iostream>

#include "barycentric_dot_product.h"
#include "global_variables.h"

/*
  this function converts a barycentric displacement vector to a cartesian vector (represented by a complex number) in the tangent space of the face

  Inputs
  u: 3d barycentric displacement vector (so u must sum to 0)
  l: |F|x3 edge lengths for each face side
  f: the face containing barycentric displacement vector u

  Outputs
  a complex number representing u in f's tangent space
*/
std::complex<double> barycentric_vector_to_cartesian(
                                                     const Eigen::Vector3d & u,
                                                     const Eigen::MatrixXd & l,
                                                     const int & f
                                                     );
#endif

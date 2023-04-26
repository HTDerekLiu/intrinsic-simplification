#ifndef BARYCENTRIC_DOT_PRODUCT
#define BARYCENTRIC_DOT_PRODUCT

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <set>
#include <iostream>

#include "global_variables.h"

/*
  this function computes the dot product of two vectors given in barcentric coordinates

  Inputs
  u: 3d barycentric displacement vector (so u must sum to 0)
  v: 3d barycentric displacement vector (so v must sum to 0)
  l: |F|x3 edge lengths for each face side
  f: the face containing barycentric displacements u and v

  Outputs
  the dot product of u with v
*/
double barycentric_dot_product(
                               const Eigen::Vector3d & u,
                               const Eigen::Vector3d & v,
                               const Eigen::MatrixXd & l,
                               const int & f
                               );
#endif

#ifndef FACE_SIDE_VECTORS_IN_FACE
#define FACE_SIDE_VECTORS_IN_FACE

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <set>
#include <iostream>

#include "opposite_corner_angle.h"
#include "global_variables.h"

/*
  this function computes the vector corresponding to a face side in its face's tangent space

  Inputs
  l:  |F|x3 edge lengths for each face side

  Outputs
  fsv: |F|x3 complex matrix of tangent vectors for each face side
*/
void face_side_vectors_in_face(
                               const Eigen::MatrixXd & l,
                               Eigen::MatrixXcd & fsv
                               );
#endif

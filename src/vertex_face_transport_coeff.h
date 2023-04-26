#ifndef VERTEX_FACE_TRANSPORT_COEFF
#define VERTEX_FACE_TRANSPORT_COEFF

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <set>
#include <iostream>

#include "global_variables.h"

/*
  this function computes the rotation taking a vertex's tangent space to the tangent space of an adjacent face

  Inputs
  fsv: |F|x3 complex matrix of tangent vectors to each face side in its face's tangent space
  A:   |F|x3 array of angular coordinates. A(f,s) gives you the angular coordinate of this face-side from its starting vertex

  Outputs
  tvf: |F|x3 complex matrix of transport coefficients. tvf(f, s) gives the rotation taking you from the tangent space of vertex F(f, s) to the tangent space of face f (Transporting from Vertex to Face)
*/
void vertex_face_transport_coeff(
                               const Eigen::MatrixXcd & fsv,
                               const Eigen::MatrixXd & A,
                               Eigen::MatrixXcd & tvf
                               );
#endif

#ifndef LAY_OUT_EDGE_DIAMOND
#define LAY_OUT_EDGE_DIAMOND

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <set>
#include <iostream>

#include "global_variables.h"

/*
  This function computes a mapping of an edge diamond in the plane.
  The diamond is always rotated so that f's first halfedge points along the positive x-axis

  Inputs
  l: |F|x3 edge lengths for each face side
  G: |F|x6 glue map
  f: the left face of the edge diamond
  s: the edge of f whose diamond we wish to lay out

  Outputs
  qi: position of the s'th vertex of face f
  qj: position of the s+1'st vertex of face f
  qk: position of the s+2'nd vertex of face f
  qm: position of the opposite vertex of opposite face
*/
void lay_out_edge_diamond(
                          const Eigen::MatrixXd & l,
                          const Eigen::MatrixXi & G,
                          const int & f,
                          const int & s,
                          std::complex<double> & qi,
                          std::complex<double> & qj,
                          std::complex<double> & qk,
                          std::complex<double> & qm
                         );
#endif

#ifndef GET_VERTEX_VECTOR_PROLONGATION
#define GET_VERTEX_VECTOR_PROLONGATION

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>

#include "vertex_tangent_space_correspondence.h"

/*
  compute the prolongation operator for vector fields at vertices

  ===Inputs
  original mesh:
  FO:     |F|x3 face list
  GO:     |F|x6 glue map
  lO:     |F|x3 edge lengths for each face side
  AO:     |F|x3 angular coordinate for each face EIGEN_STRIDE_H
  vO2fsO: |V|x2 where v2fs.row(i) returns a face side for that vertex

  coarse mesh:
  Fc:     |F|x3 face list
  Gc:     |F|x6 glue map
  lc:     |F|x3 edge lengths for each face side
  Ac:     |F|x3 angular coordinate for each face EIGEN_STRIDE_H
  vc2fsc: |V|x2 where v2fs.row(i) returns a face side for that vertex

  correspondence:
  vIdx:   a subset of coarse vertex indices, such that Vc = VO(vIdx,:)
  F2V:    |F|x{variable} list of which original vertices lie in each coarse face
  BC:     |V|x3 list of barycentric coordinates describing the locations of original vertices on the coarse mesh

  ===Outputs
  correspondence: |VO|x|Vc| complex matrix mapping coarse vector fields to vector fields on the original mesh
*/

void get_vertex_vector_prolongation(
                                    const Eigen::MatrixXi & FO,
                                    const Eigen::MatrixXi & GO,
                                    const Eigen::MatrixXd & lO,
                                    const Eigen::MatrixXd & AO,
                                    const Eigen::MatrixXi & vO2fsO,
                                    const Eigen::MatrixXi & Fc,
                                    const Eigen::MatrixXi & Gc,
                                    const Eigen::MatrixXd & lc,
                                    const Eigen::MatrixXd & Ac,
                                    const Eigen::MatrixXi & vc2fsc,
                                    const Eigen::VectorXi & vIdx,
                                    const std::vector<std::vector<int>> & F2V,
                                    const Eigen::MatrixXd & BC,
                                    Eigen::SparseMatrix<std::complex<double>> & prolongation
                                    );

#endif

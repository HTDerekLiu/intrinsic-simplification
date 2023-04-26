#ifndef VERTEX_TANGENT_SPACE_CORRESPONDENCE
#define VERTEX_TANGENT_SPACE_CORRESPONDENCE

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <set>
#include <iostream>

#include "barycentric_vector_to_cartesian.h"
#include "face_side_vectors_in_face.h"
#include "vertex_face_transport_coeff.h"
#include "vertex_one_ring_face_sides.h"
#include "lay_out_edge_diamond.h"

/*
  this function computes the mapping between vertex tangent spaces on the input mesh and the corresponding face/vertex tangent spaces on the coarse mesh

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
  VO2Vc:  |VO| dimensional vector giving the index of each original vertex in the coarse mesh (set to -1 for deleted vertices)
  F2V:    |Fc|x{variable} list of which original vertices lie in each coarse face
  VO2fc:  |VO| dimensional vector giving the face of the coarse mesh that each original vertex lies in (set to -1 for shared vertices)
  BC:     |VO|x3 list of barycentric coordinates describing the locations of original vertices on the coarse mesh

  optional:
  scale_vectors (default = false) : compute scaling component of correspondence in addition to rotation

  ===Outputs
  correspondence: |VO| dimensional list
    correspondence(v) gives v's input x axis in the coarse coordinate system, represented as a complex number
  If v is shared in both meshes, correspondence(v) maps from v's tangent space in the fine mesh to v's tangent space in the coarse mesh.
  Otherwise, correspondence(v) maps from v's tangent space in the fine mesh to the tangent space of the face containing v in the coarse mesh.
*/
void vertex_tangent_space_correspondence(
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
                                         const Eigen::VectorXi & VO2Vc,
                                         const std::vector<std::vector<int>> & F2V,
                                         const Eigen::VectorXi & VO2fc,
                                         const Eigen::MatrixXd & BC,
                                         Eigen::VectorXcd & correspondence,
                                         bool scale_vectors = false
                                        );
#endif

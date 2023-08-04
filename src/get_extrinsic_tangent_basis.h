#ifndef GET_EXTRINSIC_TANGENT_BASIS
#define GET_EXTRINSIC_TANGENT_BASIS

#include <Eigen/Core>

#include <igl/per_vertex_normals.h>

/*
  compute a tangent basis for each vertex which represents the provided angular coordinates as well as possible

Inputs:
  V: |V|x3 list of vertex positions
  F: |F|x3 vertex-face adjacency list
  A: |F|x3 angular coordinate for each face side
  v2fs: |V|x2 where v2fs.row(i) returns a face side for that vertex
Outputs:
  basisX: |V|x3 list of first tangent space basis vectors
  basisY: |V|x3 list of second tangent space basis vectors
*/

void get_extrinsic_tangent_basis(
                                 const Eigen::MatrixXd & V,
                                 const Eigen::MatrixXi & F,
                                 const Eigen::MatrixXd & A,
                                 const Eigen::MatrixXi & v2fs,
                                 Eigen::MatrixXd & basisX,
                                 Eigen::MatrixXd & basisY);
#endif

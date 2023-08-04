#ifndef GET_CONNECTION_LAPLACIAN
#define GET_CONNECTION_LAPLACIAN

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <complex>

#include <igl/cotmatrix_entries.h>

/*
  computes the vertex connection Laplacian as a complex matrix

Inputs:
  F: |F|x3 vertex-face adjacency list
  G: |F|x6 glue map
  l: |F|x3 edge lengths for each face side
Optional Inputs:
  A: |F|x3 angular coordinate for each face side
  v2fs: |v|x2 where v2fs.row(i) returns a face side for vertex i
  (if not specified, these will be computed by build_angular_coordinates)
Outputs:
  L: |V|x|V| complex connection Laplacian
*/

void get_vertex_connection_laplacian(
                                     const Eigen::MatrixXi & F,
                                     const Eigen::MatrixXi & G,
                                     const Eigen::MatrixXd & l,
                                     Eigen::SparseMatrix<std::complex<double>> & L);

void get_vertex_connection_laplacian(
                                     const Eigen::MatrixXi & F,
                                     const Eigen::MatrixXi & G,
                                     const Eigen::MatrixXd & l,
                                     const Eigen::MatrixXd & A,
                                     const Eigen::MatrixXi & v2fs,
                                     Eigen::SparseMatrix<std::complex<double>> & L);
#endif

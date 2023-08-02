#ifndef BUILD_INTRINSIC_DATA
#define BUILD_INTRINSIC_DATA

#include <Eigen/Core>
#include <Eigen/Dense>

#include "build_glue_map.h"
#include "build_edge_lengths.h"
#include "build_angular_coordinates.h"

/*
build the information needed for intrinsic triangulation 
Inputs:
    V: |V|x3 vertex list
    F: |F|x3 face list
Outputs:
    G: |F|x6 glue map
    l: |F|x3 edge lengths for each face side
    A: |F|x3 angular coordinate for each face EIGEN_STRIDE_H
    v2fs: |v|x2 where v2fs.row(i) returns a face side for that vertex
*/
void build_intrinsic_info(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs);
#endif
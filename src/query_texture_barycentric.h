#ifndef QUERY_TEXTURE_BARYCENTRIC
#define QUERY_TEXTURE_BARYCENTRIC

#include <Eigen/Core>
#include <Eigen/Dense>

#include <assert.h>
#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <tuple>

#include <igl/AABB.h>
#include "compute_barycentric_robust.h"
#include "global_variables.h"
/*
  get the barycentric coordinate for each pixel on the UV mesh (UV, UF)

Inputs:
    UV: |TC| x 2 matrix of texture coordinates (|TC| is the number of distinct texture coordinates)
    UF: |F| x 3 matrix of face indices into texture coordinate matrix
    tex_width: desired texture width in pixels (generated textures are always square)

Outputs:
    bary_faces: tex_width^2 vector indicating which face the center each pixel lies in
    bary_cords: tex_width^2 x 3 matrix of barycentric coordinates for the center of each pixel in its face
    hit_mask:   tex_width^2 vector indicating whether or not each pixel lies in a face

*/
void query_texture_barycentric(
                               const Eigen::MatrixXd & UV,
                               const Eigen::MatrixXi & UF,
                               const int & tex_width,
                               Eigen::VectorXi & bary_faces,
                               Eigen::MatrixXd & bary_coords,
                               Eigen::Matrix<bool, Eigen::Dynamic, 1> & hit_mask);

#endif

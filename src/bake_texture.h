#ifndef BAKE_TEXTURE
#define BAKE_TEXTURE

#include <Eigen/Core>
#include <Eigen/Dense>

#include <utility>
#include <array>
#include <vector>
#include <iostream>
#include <set>
#include <math.h>
#include <string>

#include <polyscope/../../deps/stb/stb_image_write.h>

/*
  These functions build and save texture maps to visualize a simplified intrinsic triangulation

  Inputs:
  png_path: file to write texture to
  F: |F|x3 vertex-face adjacency list, where |F| is the number of faces in the simplified mesh
  F2V: |F| length list of lists, where F2V[f] gives you a list of indices in BC. For example, if F2V[f] = [v], then BC[v,:] corresponds to the barycentric coordinates in F[f,:]
  hit_mask: binary mask recording whether or not each pixel lies on the parameterized mesh
  nV: number of vertices in the simplified mesh

  image_buffer: buffer to write image to, rather than saving to a file
*/

void bake_texture(
    const std::string & png_path,
    const Eigen::MatrixXi & F,
    const std::vector<std::vector<int>> & F2V,
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> & hit_mask,
    const int & nV);

void bake_texture(
    std::vector<unsigned char> & image_buffer,
    const Eigen::MatrixXi & F,
    const std::vector<std::vector<int>> & F2V,
    const Eigen::Matrix<bool, Eigen::Dynamic, 1> & hit_mask,
    const int & nV);

// image_rgba_buffer does not get modified, but the stbi_image_write function is not marked const
void bake_texture(
    const std::string & png_path,
    std::vector<unsigned char> & image_rgba_buffer);

#endif

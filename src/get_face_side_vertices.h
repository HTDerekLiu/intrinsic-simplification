#ifndef GET_FACE_SIDE_VERTICES
#define GET_FACE_SIDE_VERTICES

#include <Eigen/Core>
#include <Eigen/Dense>

// Get the vertex indices vi, vj for face f at side s
// Inputs
//     F: |F|x3 array of face list
//     fs: face-side
// Outputs
//     vi: tail vertex index of fs
//     vj: tip vertex index of fs
void get_face_side_vertices(
    const Eigen::MatrixXi & F,
    const Eigen::Vector2i & fs,
    int & vi,
    int & vj);
#endif
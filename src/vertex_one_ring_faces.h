#ifndef VERTEX_ONE_RING_FACES
#define VERTEX_ONE_RING_FACES

#include <Eigen/Core>
#include <Eigen/Dense>

#include "vertex_one_ring_unique_face_sides.h"

/*
Given a vertex v, this function outputs its one-ring face indices

Inputs
    G: |F|x6 array of gluing map. 
    v2fs: |v|x2 where v2fs.row(i) returns a face side for vertex i
    v: vertex index
Outputs
    one_ring_faces: vector of one-ring face indices
*/
void vertex_one_ring_faces(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::VectorXi & one_ring_faces);

#endif
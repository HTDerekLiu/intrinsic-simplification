#ifndef VERTEX_ONE_RING_SIZE
#define VERTEX_ONE_RING_SIZE

#include <Eigen/Core>
#include <Eigen/Dense>

#include "vertex_one_ring_face_sides.h"
#include "get_face_side_vertices.h"
#include "is_same_face_side.h"
#include "global_variables.h"
#include "ccw.h"

/*
Given a vertex v, this function outputs the number of its one-ring vertices

Inputs
    G: |F|x6 array of gluing map. 
    v2fs: |v|x2 where v2fs.row(i) returns a face side for vertex i
    v: vertex index
Outputs
    size of the one-ring
*/
int vertex_one_ring_size(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v);

#endif
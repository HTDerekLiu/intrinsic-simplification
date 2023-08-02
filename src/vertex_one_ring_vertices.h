#ifndef VERTEX_ONE_RING_VERTICES
#define VERTEX_ONE_RING_VERTICES

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <unordered_set>

#include "vertex_one_ring_face_sides.h"
#include "get_face_side_vertices.h"
#include "is_same_face_side.h"
#include "global_variables.h"
#include "ccw.h"

/*
Given a vertex v, this function outputs its one-ring vertex indices

Inputs
    G: |F|x6 array of gluing map. 
    G: |F|x6 glue map
    v2fs: |v|x2 where v2fs.row(i) returns a face side for vertex i
    v: vertex index
Outputs
    one_ring_vertices: vector of one-ring vertex indices
    is_boundary_vertex: bool variable to show whether v is the boundary vertex
    face_sides_of_one_ring_vertices: as the name stated

*/
void vertex_one_ring_vertices(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::VectorXi & one_ring_vertices,
    bool & is_boundary_vertex,
    std::vector<Eigen::Vector2i> & face_sides_of_one_ring_vertices);

void vertex_one_ring_vertices(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::VectorXi & one_ring_vertices,
    bool & is_boundary_vertex);

void vertex_one_ring_vertices(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::VectorXi & one_ring_vertices);

#endif
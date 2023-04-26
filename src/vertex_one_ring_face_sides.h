#ifndef VERTEX_ONE_RING_FACE_SIDES
#define VERTEX_ONE_RING_FACE_SIDES

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include "next.h"
#include "twin.h"
#include "ccw.h"
#include "cw.h"
#include "is_same_face_side.h"
#include "global_variables.h"
#include "get_face_side_vertices.h"

/*
    Get the vertex one-ring face-sides starting from the vertes
    Inputs
        G: |F|x6 array of gluing map. 
        fs_start: starting face side
    Outputs
        one_ring_fs: std::vector of one-ring face sides

    Note
    if the vertex v is a boundary vertex, this will only return the valid face-sides in the ccw order. Therefore:
    - one_ring_fs.at(0).twin is GHOST_FACE_SIDE
    - one_ring_fs.at(-1).next.next.twin is also a GHOST_FACE_SIDE. 

    [WARNING]
    This function will return all the face sides of a vertex, including two half-edges of a self-edge. Therefore, one_ring_fs.size() != the number of actual one-ring faces. 
*/

void vertex_one_ring_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs_start,
    std::vector<Eigen::Vector2i> & one_ring_fs);

// Get the vertex one-ring face-sides starting from the vertes
// Inputs
//     G: |F|x6 array of gluing map. 
//     fs_start: starting face side
// Outputs
//     one_ring_fs: std::vector of one-ring face sides
//     is_boundary_vertex: boolean variable indicating whether the vertex at the tail of this fs is on the boundary 
void vertex_one_ring_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs_start,
    std::vector<Eigen::Vector2i> & one_ring_fs,
    bool & is_boundary_vertex);

void vertex_one_ring_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs, // vertex to fs map
    const int & v, // the vertex you want to extract one-ring
    std::vector<Eigen::Vector2i> & one_ring_fs);

void vertex_one_ring_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs, // vertex to fs map
    const int & v, // the vertex you want to extract one-ring
    std::vector<Eigen::Vector2i> & one_ring_fs,
    bool & is_boundary_vertex);
#endif
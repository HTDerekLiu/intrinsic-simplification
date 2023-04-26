#ifndef VERTEX_ONE_RING_UNIQUE_FACE_SIDES
#define VERTEX_ONE_RING_UNIQUE_FACE_SIDES

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <unordered_set>
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
    one_ring_fs: std::vector of one-ring face sides with unique faces

Note #1
if the vertex v is a boundary vertex, this will only return the valid face-sides in the ccw order. Therefore:
- one_ring_fs.at(0).twin is GHOST_FACE_SIDE
- one_ring_fs.at(-1).next.next.twin is also a GHOST_FACE_SIDE. 

Note #2
compared to "vertex_one_ring_face_sides.h" this function DOES NOT return all the face sides from this vertex. Instead, it only returns a list of face sides with unique faces. In 99% cases, "vertex_one_ring_face_sides" and "vertex_one_ring_unique_face_sides" will give you the same result, the only situation it will be different is the existance of self-edge (an edge connecting to the same vertex)
*/
void vertex_one_ring_unique_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs_start,
    std::vector<Eigen::Vector2i> & one_ring_fs);

// Get the vertex one-ring face-sides starting from the vertes
// Inputs
//     G: |F|x6 array of gluing map. 
//     fs_start: starting face side
// Outputs
//     one_ring_fs: std::vector of one-ring face sides with unique faces
//     is_boundary_vertex: boolean variable indicating whether the vertex at the tail of this fs is on the boundary 
void vertex_one_ring_unique_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::Vector2i & fs_start,
    std::vector<Eigen::Vector2i> & one_ring_fs,
    bool & is_boundary_vertex);

void vertex_one_ring_unique_face_sides(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs, // vertex to fs map
    const int & v, // the vertex you want to extract one-ring
    std::vector<Eigen::Vector2i> & one_ring_fs);
#endif
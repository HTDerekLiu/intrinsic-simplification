#ifndef GAUSSIAN_CURVATURE_AT_VERTEX
#define GAUSSIAN_CURVATURE_AT_VERTEX

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

#include <opposite_corner_angle.h>
#include <twin.h>
#include <next.h>
#include <is_same_face_side.h>
#include <is_interior_vertex.h>
#include <vertex_one_ring_unique_face_sides.h>
#include <vertex_one_ring_face_sides.h>
#include <global_variables.h>
/*
Computes discrete gaussian curvature for a given vertex. For boundary vertices, we output geodesic curvature

Input
G: |F|x3x2 array of gluing map. 
El: |F|x3 array of face-side lengths
v2fs: |V|x2 array of face-sides. V2FS[v] outputs one face-side starting from this vertex 
v: vertex index

Outputs
gaussian curvature at v
*/
double gaussian_curvature_at_vertex(
    const Eigen::MatrixXi &G, 
    const Eigen::MatrixXd &l, 
    const Eigen::MatrixXi &v2fs, 
    const int vertex);

#endif
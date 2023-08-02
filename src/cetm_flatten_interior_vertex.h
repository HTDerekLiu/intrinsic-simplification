#ifndef CETM_FLATTEN_INTERIOR_VERTEX
#define CETM_FLATTEN_INTERIOR_VERTEX

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <math.h> 
#include <set>

#include "is_interior_vertex.h"
#include "vertex_one_ring_face_sides.h"
#include "opposite_corner_angle.h"
#include "vertex_one_ring_vertices.h"
#include "vertex_one_ring_faces.h"
#include "get_face_side_vertices.h"
#include "gaussian_curvature_at_vertex.h"
#include "polar_to_cartesian.h"

// This is a helper function for "flatten_interior_vertex_and_cost.h", which flatten an interior vertex with CETM. One should use "flatten_interior_vertex_and_cost" instead
void cetm_flatten_interior_vertex(
    const Eigen::MatrixXd & l012,
    const Eigen::MatrixXd & s012,
    double & u);

#endif
#ifndef COARSEN_MESH
#define COARSEN_MESH

#include <Eigen/Core>
#include <Eigen/Dense>

#include <assert.h>
#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <tuple>

#include "global_variables.h"
#include "gaussian_curvature_at_vertex.h"
#include "flatten_interior_vertex_and_cost.h"
#include "flatten_boundary_vertex_and_cost.h"
#include "flatten_ear_vertex_and_cost.h"
#include "always_flip_to_degree_three.h"
#include "always_flip_to_ear.h"
#include "remove_degree_three_vertex.h"
#include "remove_ear_vertex.h"
#include "vertex_areas.h"
#include "is_interior_vertex.h"
#include "is_ear_vertex.h"
#include "is_glue_map_valid.h"

/*
intrinsic coarsening
*/

void coarsen_mesh(
    const int & total_removal,
    const double & mixture_weight, // 0 is pure curvature, 1 is pure area
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V);

#endif
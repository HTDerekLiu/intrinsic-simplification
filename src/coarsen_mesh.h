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
This function simplifies a mesh intrinsically.

Inputs:
    total_removal: scalar. number of vertices one wants to remove
    mixture_weight: scalar between 0 and 1. This parameter determines which geometric quantities to use (0 for pure curvature, 1 for pure area, 0~1 uses a mixture of them)
    F: |F|x3 vertex-face adjacency list
    G: |F|x6 glue map
    l: |F|x3 edge lengths for each face side
    A: |F|x3 angular coordinate for each face side
    v2fs: |v|x2 where v2fs.row(i) returns a face side for vertex i
    BC: |BC|x3 array of barycentric coordinates whose corresponding faces are stored in F2V implicitly
    F2V: |F| length list of lists, where F2V[f] gives you a list of indices in BC. For example, if F2V[f] = [v], then BC[v,:] corresponds to the barycentric coordinates in F[f,:]
Outputs:
    F,G,l,A,v2fs,BC,F2V are changed in place
*/
void coarsen_mesh(
    const int & total_removal,
    const double & mixture_weight, 
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V);

#endif
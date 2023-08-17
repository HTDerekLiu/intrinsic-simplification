#ifndef DELAUNAY_REFINE
#define DELAUNAY_REFINE

#include <Eigen/Core>
#include <Eigen/Dense>

#include "global_variables.h"

/*
    This function performs intrinsic delaunay refinement and transfers barycentric coordinates to the refined mesh

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
    angle_threshold_degrees: minimum allowed angle in refined triangulation, specified in degrees. Values > 30 degrees may not terminate
    circumradius_threshold: maximum allowed triangle circumradius in refined triangulation
    max_iterations: maximum number of iterations in refinement algorithm. If negative, algorithm runs until quality criteria are satisfied.
Outputs:
    F, G, l, A, v2fs, BC, F2V are changed in place
*/
void delaunay_refine(
  Eigen::MatrixXi & F,
  Eigen::MatrixXi & G,
  Eigen::MatrixXd & l,
  Eigen::MatrixXd & A,
  Eigen::MatrixXi & v2fs,
  Eigen::MatrixXd & BC,
  std::vector<std::vector<int>> & F2V,
  double angle_threshold_degrees = 25,
  double circumradius_threshold = std::numeric_limits<double>::infinity(),
  int max_iterations = -1);

#endif

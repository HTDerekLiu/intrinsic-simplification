#ifndef CONNECTED_COMPONENTS
#define CONNECTED_COMPONENTS

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <queue>

#include <global_variables.h>
#include <next.h>
#include <twin.h>

/*
  CONNECTED_COMPONENTS finds the connected components of a mesh

Inputs:
  F: |F|x3 vertex-face adjacency list F
  G: |F|x6 array of gluing map.

Outputs:
  n_components: the number of components
  v_id: |V| vector of vertex component ids
  f_id: |F| vector of face component ids
*/
void connected_components(
                          const Eigen::MatrixXi & F,
                          const Eigen::MatrixXi & G,
                          int & n_components,
                          Eigen::VectorXi & v_id,
                          Eigen::VectorXi & f_id);
#endif

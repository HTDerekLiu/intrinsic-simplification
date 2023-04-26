#ifndef HEAT_GEODESICS
#define HEAT_GEODESICS

#include <igl/min_quad_with_fixed.h>
#include <Eigen/Sparse>
#include <Eigen/Sparse>

#include <igl/grad.h>
#include <igl/grad_intrinsic.h>
#include <igl/doublearea.h>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_intrinsic.h>
#include <igl/intrinsic_delaunay_cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/massmatrix_intrinsic.h>
#include <igl/grad_intrinsic.h>
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/slice.h>
#include <igl/avg_edge_length.h>

#include "build_edge_lengths.h"

// This code is adapted from libigl https://github.com/libigl/libigl/blob/main/include/igl/heat_geodesics.h
struct heat_geodesics_data
{
  // Gradient and Divergence operators
  Eigen::SparseMatrix<double> Grad,Div;
  // Number of gradient components
  int ng;
  // List of boundary vertex indices
  Eigen::VectorXi b;
  // Solvers for Dirichet, Neumann problems
  igl::min_quad_with_fixed_data<double> Dirichlet, Neumann, Poisson;
};

// Precompute factorized solvers for computing a fast approximation of
// geodesic distances on a mesh (V,F). [Crane et al. 2013]
//
// Inputs:
//   l  #F by 3 list of face side length
//   F  #F by 3 list of mesh face indices into V
// Outputs:
//   data  precomputation data (see heat_geodesics_solve)
bool heat_geodesics_precompute(
    const Eigen::MatrixXd & Vorl,
    const Eigen::MatrixXi & F,
    heat_geodesics_data & data);

// Compute fast approximate geodesic distances using precomputed data from a
// set of selected source vertices (gamma)
//
// Inputs: 
//   data  precomputation data (see heat_geodesics_precompute)
//   gamma  #gamma list of indices into V of source vertices
// Outputs:
//   D  #V list of distances to gamma 
void heat_geodesics_solve(
    const heat_geodesics_data & data,
    const Eigen::VectorXi & gamma,
    Eigen::VectorXd & D);

#endif
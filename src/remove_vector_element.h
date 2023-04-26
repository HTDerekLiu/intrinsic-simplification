#ifndef REMOVE_VECTOR_ELEMENT_H
#define REMOVE_VECTOR_ELEMENT_H

#include <Eigen/Core>

// remove a element in a Eigen VectorXd or VectorXi
// Inputs: 
//   idx: index to remove
//   vec: Eigen Vector 
// Outputs:
//   vec removed Eigen Vector
  
void remove_vector_element(
  const int idx,
  Eigen::VectorXd & vec);

void remove_vector_element(
  const int idx,
  Eigen::VectorXi & vec);

#endif
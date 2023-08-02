#ifndef SORT_VEC
#define SORT_VEC

#include<algorithm>
#include<Eigen/Core>

/*
sort vector

Inputs:
  vec: input vector

Outputs:
  sorted_vec: sorted vector
  ind: indices of the sorting (argsort)
*/

template <typename Derived>
void sort_vec(
  const Eigen::PlainObjectBase<Derived> & vec, 
  Eigen::PlainObjectBase<Derived> & sorted_vec,  
  Eigen::VectorXi & ind);
#endif
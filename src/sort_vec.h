#ifndef SORT_VEC
#define SORT_VEC

#include<algorithm>
#include<Eigen/Core>

template <typename Derived>
void sort_vec(
  const Eigen::PlainObjectBase<Derived> & vec, 
  Eigen::PlainObjectBase<Derived> & sorted_vec,  
  Eigen::VectorXi & ind);
#endif
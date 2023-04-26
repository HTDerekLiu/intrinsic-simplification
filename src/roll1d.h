#ifndef ROLL1D
#define ROLL1D

#include <Eigen/Core>
#include <Eigen/Dense>

// roll a Eigen Vector (the same as numpy.roll)
template <typename Derived>
void roll1d(
    const Eigen::PlainObjectBase<Derived> & vin,
    const int & n, 
    Eigen::PlainObjectBase<Derived> & vout);
#endif
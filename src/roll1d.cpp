#include "roll1d.h"

template <typename Derived>
void roll1d(
    const Eigen::PlainObjectBase<Derived> & vin,
    const int & n, 
    Eigen::PlainObjectBase<Derived> & vout)
{
    int length = vin.size();
    assert(n < length); // n must not bigger than the size of hte vector

    vout.tail(length-n) = vin.head(length - n);
    vout.head(n) = vin.tail(n);
}

template void roll1d<Eigen::Matrix<int, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, 3, 1, 0, 3, 1> > const&, int const&, Eigen::PlainObjectBase<Eigen::Matrix<int, 3, 1, 0, 3, 1> >&);

template void roll1d<Eigen::Matrix<double, 3, 1, 0, 3, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> > const&, int const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 3, 1, 0, 3, 1> >&);
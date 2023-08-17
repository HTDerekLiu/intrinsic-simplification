#include "face_circumradius.h"

double face_circumradius(
                         const Eigen::MatrixXd & l,
                         int iF) {
  double a = l(iF, 0);
  double b = l(iF, 1);
  double c = l(iF, 2);

  // compute area by Heron's rule
  double s = (a + b + c) / 2.0;
  double squared_area = s * (s - a) * (s - b) * (s - c);
  squared_area = std::fmax(0., squared_area); // clamp squared area to be positive
  double area = std::sqrt(squared_area);

  return a * b * c / (4. * area);
}

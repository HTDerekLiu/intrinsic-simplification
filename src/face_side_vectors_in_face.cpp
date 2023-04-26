#include "face_side_vectors_in_face.h"

void face_side_vectors_in_face(
                               const Eigen::MatrixXd & l,
                               Eigen::MatrixXcd & fsv
                               )
{
  fsv.resize(l.rows(), 3);
  for (int f=0; f<l.rows(); f++)
  {
    double l0 = l(f, 0);
    double l1 = l(f, 1);
    double l2 = l(f, 2);

    double theta1 = opposite_corner_angle(l, Eigen::Vector2i{f, 1});
    double theta2 = opposite_corner_angle(l, Eigen::Vector2i{f, 2});

    fsv(f, 0) = std::polar(l0, 0.);
    fsv(f, 1) = std::polar(l1, M_PI - theta2);
    fsv(f, 2) = std::polar(l2, theta1 - M_PI);
  }
}

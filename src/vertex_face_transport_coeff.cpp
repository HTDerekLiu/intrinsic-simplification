#include "vertex_face_transport_coeff.h"

void vertex_face_transport_coeff(
                                 const Eigen::MatrixXcd & fsv,
                                 const Eigen::MatrixXd & A,
                                 Eigen::MatrixXcd & tvf
                                 )
{
  tvf.resize(fsv.rows(), 3);
  for (int f=0; f<fsv.rows(); f++)
  {
    for (int s=0; s<3; s++)
    {
      std::complex<double> fs_vertex = std::polar(1., A(f, s));
      std::complex<double> fs_face = fsv(f, s) / std::abs(fsv(f, s));

      tvf(f, s) = fs_face / fs_vertex;
    }
  }
}

#include "lay_out_edge_diamond.h"

void lay_out_edge_diamond(
                          const Eigen::MatrixXd & l,
                          const Eigen::MatrixXi & G,
                          const int & f,
                          const int & s,
                          std::complex<double> & qi,
                          std::complex<double> & qj,
                          std::complex<double> & qk,
                          std::complex<double> & qm
                         )
{

  auto triangle_area = [&](double lij, double ljk, double lki) {
                         double s = (lij + ljk + lki) / 2.;
                         return sqrt(abs(s*(s-lij)*(s-ljk)*(s-lki)));
                       };
  int g = G(f, 2*s);
  int sg = G(f, 2*s+1);

  // face f is ijk, g is jim

  double lij = l(f, (s+0)%3);
  double ljk = l(f, (s+1)%3);
  double lki = l(f, (s+2)%3);
  double Aijk = triangle_area(lij, ljk, lki);

  double ijk_acute = (ljk * ljk < lij * lij + lki * lki) ? 1 : -1;
  double ky = 2. * Aijk / lij;
  double kx = ijk_acute * sqrt(abs(lki * lki - ky * ky));

  double lim = l(g, (sg+1)%3);
  double lmj = l(g, (sg+2)%3);
  double Ajim = triangle_area(lij, lim, lmj);

  bool imj_acute = (lmj * lmj < lij * lij + lim * lim) ? 1 : -1;
  double my = -2. * Ajim / lij;
  double mx = imj_acute * sqrt(abs(lim * lim - my * my));

  qi = std::complex<double>{0., 0.};
  qj = std::complex<double>{lij, 0.};
  qk = std::complex<double>{kx, ky};
  qm = std::complex<double>{mx, my};

  std::complex<double> X_ijk; // tangent basis vector for face f
  if (s==0) {
    X_ijk = (qj - qi);
  } else if (s == 1) {
    X_ijk = (qi - qk);
  } else if (s == 2) {
    X_ijk = (qk - qj);
  }
  X_ijk /= std::abs(X_ijk);

  qi /= X_ijk;
  qj /= X_ijk;
  qk /= X_ijk;
  qm /= X_ijk;
}

#include "get_vertex_vector_prolongation.h"

void get_vertex_vector_prolongation(
                                    const Eigen::MatrixXi & FO,
                                    const Eigen::MatrixXi & GO,
                                    const Eigen::MatrixXd & lO,
                                    const Eigen::MatrixXd & AO,
                                    const Eigen::MatrixXi & vO2fsO,
                                    const Eigen::MatrixXi & Fc,
                                    const Eigen::MatrixXi & Gc,
                                    const Eigen::MatrixXd & lc,
                                    const Eigen::MatrixXd & Ac,
                                    const Eigen::MatrixXi & vc2fsc,
                                    const Eigen::VectorXi & vIdx,
                                    const std::vector<std::vector<int>> & F2V,
                                    const Eigen::MatrixXd & BC,
                                    Eigen::SparseMatrix<std::complex<double>> & prolongation
                                    )
{
  using namespace Eigen;

  // ==== init
  VectorXi VO2fc = VectorXi::Constant(vO2fsO.rows(), -1);
  for (int fc=0; fc<F2V.size(); fc++)
    {
      for (int vO : F2V[fc])
        {
          if (vO > VO2fc.rows()) continue;
          VO2fc(vO) = fc;
        }
    }

  VectorXi VO2Vc = VectorXi::Constant(vO2fsO.rows(), -1);
  for (int vc=0; vc<vc2fsc.rows(); vc++)
    VO2Vc(vIdx(vc)) = vc;


  VectorXcd correspondence;
  vertex_tangent_space_correspondence(FO, GO, lO, AO, vO2fsO, Fc, Gc, lc, Ac, vc2fsc, vIdx, VO2Vc, F2V, VO2fc, BC, correspondence);

  // auto toCplx = [&](VectorXd v)->std::complex<double>{return std::complex<double>{v(0), v(1)};};

  // auto vectorBarycentricInterpolationMap = [&](int f, VectorXd bc) -> std::complex<double> {
  //   std::complex<double> xi_i = toCplx(X.row(F(f, 0)));
  //   std::complex<double> xj_j = toCplx(X.row(F(f, 1)));
  //   std::complex<double> xk_k = toCplx(X.row(F(f, 2)));

  //   std::complex<double> xi_f = tvf(f, 0) * xi_i;
  //   std::complex<double> xj_f = tvf(f, 1) * xj_j;
  //   std::complex<double> xk_f = tvf(f, 2) * xk_k;

  //   std::complex<double> x = bc(0) * xi_f + bc(1) * xj_f + bc(2) * xk_f;
  //   return x;
  // };

  MatrixXcd fsv;
  face_side_vectors_in_face(lc, fsv);

  MatrixXcd tvf; // TODO: symmetrize transport coeff
  vertex_face_transport_coeff(fsv, Ac, tvf);

  //== define coarse vector field and interpolate across faces
  std::vector<Triplet<std::complex<double>>> T;

  // Prolong
  for (int vO=0; vO<vO2fsO.rows(); vO++){
    int fc = VO2fc(vO);
    int vc = VO2Vc(vO);

    // HACK: skip prolonging vectors if we failed to find correspondence
    if (std::abs(correspondence(vO)) < 1e-8) continue;

    std::complex<double> x_coarse;
    if (vc >= 0) { // shared vertex
      T.emplace_back(vO, vc, 1./correspondence(vO));
    } else {
      std::complex<double> c0 = BC(vO, 0) * tvf(fc, 0) / correspondence(vO);
      std::complex<double> c1 = BC(vO, 1) * tvf(fc, 1) / correspondence(vO);
      std::complex<double> c2 = BC(vO, 2) * tvf(fc, 2) / correspondence(vO);

      T.emplace_back(vO, Fc(fc, 0), c0);
      T.emplace_back(vO, Fc(fc, 1), c1);
      T.emplace_back(vO, Fc(fc, 2), c2);
    }
  }

  prolongation.resize(vO2fsO.rows(), vc2fsc.rows());
  prolongation.setFromTriplets(T.begin(), T.end());
}

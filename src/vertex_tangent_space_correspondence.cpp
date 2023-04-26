#include "vertex_tangent_space_correspondence.h"

void vertex_tangent_space_correspondence(
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
                                         const Eigen::VectorXi & VO2Vc,
                                         const std::vector<std::vector<int>> & F2V,
                                         const Eigen::VectorXi & VO2fc,
                                         const Eigen::MatrixXd & BC,
                                         Eigen::VectorXcd & correspondence,
                                         bool scale_vectors
                                        )
{
  using namespace Eigen;

  correspondence.resize(vO2fsO.rows());

  auto tipVertexO = [&](VectorXi fs) -> int {return FO(fs(0), (fs(1)+1)%3);};

  MatrixXcd fsv;
  face_side_vectors_in_face(lc, fsv);

  MatrixXcd tvf; // TODO: symmetries transport coeff
  vertex_face_transport_coeff(fsv, Ac, tvf);

  // find correspondence
  for (int vO=0; vO<vO2fsO.rows(); vO++){
    int fc_v = VO2fc(vO);

    if (fc_v < 0) // vO is coarse, just use its angular coordinate
    {
      // TODO: estimate scale?
      correspondence(vO) = std::complex<double>{1., 0.};
      continue;
    }

    // if vO is not coarse, look at its neighbors
    std::vector<Vector2i> vOutgoingFaceSides;
    vertex_one_ring_face_sides(GO, vO2fsO, vO, vOutgoingFaceSides);

    std::complex<double> totalSignpost{0,0};
    double goodDegree = 0;
    for (auto fsO : vOutgoingFaceSides){
      int wO = tipVertexO(fsO);

      //== angle in fine mesh coordinates
      double theta_vwO = AO(fsO(0), fsO(1));
      double lvwO = lO(fsO(0), fsO(1));
      std::complex<double> vwO = std::polar(lvwO, theta_vwO);

      int fc_w = VO2fc(wO);

      if (fc_w == fc_v){
        // v and w lie in the same coarse face: connect them directly

        //== angle in coarse mesh coordinates
        std::complex<double> vwc = barycentric_vector_to_cartesian(BC.row(wO)-BC.row(vO), lc, fc_w);

        std::complex<double> vw_signpost = vwc / vwO;
        totalSignpost += vw_signpost;
        goodDegree += 1;
      } else {
        // check if fc_w neighbors fc_v
        for (int sv=0; sv<3; sv++) {
          if (Gc(fc_v, 2*sv) == fc_w) {
            int sw = Gc(fc_v, 2*sv+1);
            std::complex<double> qi, qj, qk, qm;
            lay_out_edge_diamond(lc, Gc, fc_v, sv, qi, qj, qk, qm);

            VectorXd vBC = BC.row(vO);
            VectorXd wBC = BC.row(wO);
            std::complex<double> qv = vBC((sv+0)%3) * qi + vBC((sv+1)%3) * qj + vBC((sv+2)%3) * qk;
            std::complex<double> qw = wBC((sw+0)%3) * qj + wBC((sw+1)%3) * qi + wBC((sw+2)%3) * qm;

            std::complex<double> vwc = (qw - qv) ;

            std::complex<double> vw_signpost = vwc / vwO;
            totalSignpost += vw_signpost;
            goodDegree += 1;
          }
        }
      }
    }

    if (goodDegree > 0){
      if (scale_vectors) {
        // preserve scaling component of transformation
        totalSignpost /= goodDegree;
      } else {
        // return just the rotation
        totalSignpost /= std::abs(totalSignpost);
      }
      correspondence(vO) = totalSignpost;
    } else {
      correspondence(vO) = std::complex<double>{0.,0.};
    }
  }
}

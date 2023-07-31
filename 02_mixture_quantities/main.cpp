#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <igl/read_triangle_mesh.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <build_intrinsic_info.h>
#include <chrono>
#include <coarsen_mesh.h>
#include <get_barycentric_points.h>
#include <remove_unreferenced_intrinsic.h>

int main(int argc, char* argv[]) {
  using namespace Eigen;
  using namespace std;
  using namespace global_variables;
  using namespace std::chrono;

  // load mesh
  MatrixXd VO;
  MatrixXi FO;
  {
    igl::read_triangle_mesh("../../meshes/capsule.obj", VO, FO);
  }
  
  int total_removal = 10000;

  // pure curvature
  MatrixXd V_curvature;
  MatrixXi F_curvature = FO;
  {
    MatrixXi G;    // glue map
    MatrixXd l;    // edge lengths
    MatrixXd A;    // angular coordinates
    MatrixXi v2fs; // vertex to faceside map
    build_intrinsic_info(VO, FO, G, l, A, v2fs);

    double weight = 0.0; // pure curvature
    MatrixXd BC;
    vector<vector<int>> F2V;
    coarsen_mesh(total_removal, weight, F_curvature, G, l, A, v2fs, BC, F2V);

    map<int, int> IMV, IMF;
    VectorXi vIdx, fIdx;
    remove_unreferenced_intrinsic(F_curvature, IMV, IMF, vIdx, fIdx);
    igl::slice(VO,vIdx,1,V_curvature);
  }

  // mixture of curvature and area
  MatrixXd V_mix;
  MatrixXi F_mix = FO;
  {
    MatrixXi G;    // glue map
    MatrixXd l;    // edge lengths
    MatrixXd A;    // angular coordinates
    MatrixXi v2fs; // vertex to faceside map
    build_intrinsic_info(VO, FO, G, l, A, v2fs);

    double weight = 0.8; // part curvature + part area. Note that they are in different units so the percentage doesn't have much physical meaning. 
    MatrixXd BC;
    vector<vector<int>> F2V;
    coarsen_mesh(total_removal, weight, F_mix, G, l, A, v2fs, BC, F2V);

    map<int, int> IMV, IMF;
    VectorXi vIdx, fIdx;
    remove_unreferenced_intrinsic(F_mix, IMV, IMF, vIdx, fIdx);
    igl::slice(VO,vIdx,1,V_mix);
  }

  // full area
  MatrixXd V_area;
  MatrixXi F_area = FO;
  {
    MatrixXi G;    // glue map
    MatrixXd l;    // edge lengths
    MatrixXd A;    // angular coordinates
    MatrixXi v2fs; // vertex to faceside map
    build_intrinsic_info(VO, FO, G, l, A, v2fs);

    double weight = 1.0; // pure area
    MatrixXd BC;
    vector<vector<int>> F2V;
    coarsen_mesh(total_removal, weight, F_area, G, l, A, v2fs, BC, F2V);

    map<int, int> IMV, IMF;
    VectorXi vIdx, fIdx;
    remove_unreferenced_intrinsic(F_area, IMV, IMF, vIdx, fIdx);
    igl::slice(VO,vIdx,1,V_area);
  }


  polyscope::init();
  polyscope::registerSurfaceMesh("input mesh", VO, FO);
  polyscope::registerSurfaceMesh("curvature metric (with wrong edge length)", V_curvature, F_curvature);
  polyscope::registerSurfaceMesh("mixture metric (with wrong edge length)", V_mix, F_mix);
  polyscope::registerSurfaceMesh("area metric (with wrong edge length)", V_area, F_area);
  polyscope::view::lookAt(glm::vec3{1.5, 1.5, 3}, glm::vec3{0., 0., 0.});
  polyscope::show();
}

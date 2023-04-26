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
  MatrixXi F, FO;
  {
    // igl::read_triangle_mesh("../../meshes/spot_subdiv.obj", VO, FO);
    igl::read_triangle_mesh("../../meshes/bs_rest.obj", VO, FO);
  }

  MatrixXi G;    // glue map
  MatrixXd l;    // edge lengths
  MatrixXd A;    // angular coordinates
  MatrixXi v2fs; // vertex to faceside map
  build_intrinsic_info(VO, FO, G, l, A, v2fs);
  F = FO;

  auto start = high_resolution_clock::now();

  int total_removal = 19500;
  double weight = 0.0; // 0: pure curvature, 1: pure area
  MatrixXd BC;
  vector<vector<int>> F2V;
  coarsen_mesh(total_removal, weight, F, G, l, A, v2fs, BC, F2V);

  auto stop     = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);

  cout << "removed " << total_removal << " vertices, it takes " << ((float) duration.count() / 1000000.0) << " seconds" << endl;

  MatrixXd P;
  get_barycentric_points(VO, F, BC, F2V, P);

  // removed unreferenced vertices
  map<int, int> IMV, IMF;
  VectorXi vIdx, fIdx;
  remove_unreferenced_intrinsic(F,IMV, IMF, vIdx, fIdx);
  MatrixXd V;
  igl::slice(VO,vIdx,1,V);

  polyscope::init();
  polyscope::registerSurfaceMesh("input mesh", VO, FO);
  polyscope::registerSurfaceMesh("coarsened mesh (with wrong edge length)", V, F);
  polyscope::registerPointCloud("barycentric points from the input", P);
  polyscope::view::lookAt(glm::vec3{1.5, 1.5, 3}, glm::vec3{0., 0., 0.});
  polyscope::show();
}

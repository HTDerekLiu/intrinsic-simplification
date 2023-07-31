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
#include <get_prolongation.h>

int main(int argc, char* argv[]) {
  using namespace Eigen;
  using namespace std;
  using namespace global_variables;
  using namespace std::chrono;

  // load mesh
  MatrixXd VO;
  MatrixXi F, FO;
  {
    igl::read_triangle_mesh("../../meshes/spot.obj", VO, FO);
  }

  MatrixXi G;    // glue map
  MatrixXd l;    // edge lengths
  MatrixXd A;    // angular coordinates
  MatrixXi v2fs; // vertex to faceside map
  build_intrinsic_info(VO, FO, G, l, A, v2fs);
  F = FO;

  int total_removal = 2000;
  double weight = 0.0; // 0: pure curvature, 1: pure area
  MatrixXd BC;
  vector<vector<int>> F2V;
  coarsen_mesh(total_removal, weight, F, G, l, A, v2fs, BC, F2V);

  cout << "removed " << total_removal << " vertices" << endl;

  // removed unreferenced vertices
  map<int, int> IMV, IMF;
  VectorXi vIdx, fIdx;
  remove_unreferenced_intrinsic(F, G, l, A, v2fs, F2V, IMV, IMF, vIdx, fIdx);
  MatrixXd V;
  igl::slice(VO,vIdx,1,V);

  // get prolongation
	SparseMatrix<double> P;
	get_prolongation(F,BC,F2V,IMV,P);

	// create a random test function
	VectorXd D(V.rows());
  D.setZero();
  D(0) = 1.0;

  // prolongation
	VectorXd DO = P * D;

  polyscope::init();
  polyscope::registerSurfaceMesh("input mesh", VO, FO);
  polyscope::getSurfaceMesh("input mesh")->addVertexScalarQuantity("prolonged rand function", DO);
  polyscope::registerSurfaceMesh("coarsened mesh (with wrong edge length)", V, F);
  polyscope::getSurfaceMesh("coarsened mesh (with wrong edge length)")->addVertexScalarQuantity("rand function", D);
  polyscope::view::lookAt(glm::vec3{1.5, 1.5, 3}, glm::vec3{0., 0., 0.});
  polyscope::show();
}

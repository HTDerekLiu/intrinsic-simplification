#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "args/args.hxx"

#include <igl/read_triangle_mesh.h>
#include <igl/cotmatrix_intrinsic.h>
#include <igl/massmatrix_intrinsic.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <build_intrinsic_info.h>
#include <chrono>
#include <coarsen_mesh.h>
#include <get_barycentric_points.h>
#include <remove_unreferenced_intrinsic.h>
#include <get_prolongation.h>

// Stolen from https://github.com/nmwsharp/nonmanifold-Laplacian src/main.cpp
void saveMatrix(Eigen::SparseMatrix<double>& matrix, std::ostream& out) {
  out << std::setprecision(16);

  for (int k = 0; k < matrix.outerSize(); ++k) {
    for (typename Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it;
         ++it) {
      double val       = it.value();
      size_t iRow = it.row();
      size_t iCol = it.col();

      out << (iRow + 1) << " " << (iCol + 1) << " " << val << std::endl;
    }
  }
}

void saveMatrix(Eigen::SparseMatrix<double>& matrix, std::string filename) {
  // WARNING: this follows matlab convention and thus is 1-indexed

  std::cout << "Writing sparse matrix to: " << filename << std::endl;

  std::ofstream outFile(filename);
  if (!outFile) {
    throw std::runtime_error("failed to open output file " + filename);
  }
  saveMatrix(matrix, outFile);
  outFile.close();
}

int main(int argc, char* argv[]) {
  using namespace Eigen;
  using namespace std;
  using namespace global_variables;
  using namespace std::chrono;

  // Configure the argument parser
  args::ArgumentParser parser("Intrinsic Prolongation");
  args::Positional<std::string> filename_arg(parser, "mesh",
                                             "Mesh to be coarsened. (default='../../meshes/spot.obj')");
  args::Positional<int> n_coarse_vertices_arg(parser, "n_coarse_vertices",
                                              "Number of coarse vertices to leave. (default='500')");
  args::ValueFlag<double> weight_arg(parser, "area_weight",
                                     "Influence of vertex area on coarsening. 0: none, 1: pure area weighting. (default='0')", {'w',"area_weight"});
  args::ValueFlag<std::string> prolongation_matrix_path_arg(parser, "prolongation_matrix_path",
                                     "File to save prolongation matrix to. If not set, the prolongation matrix is not saved", {'p',"prolongation_path"});
  args::ValueFlag<std::string> laplace_matrix_path_arg(parser, "laplace_matrix_path",
                                                            "File to save coarsened laplace matrix to. If not set, the laplace matrix is not computed", {'l',"laplace_path"});
  args::ValueFlag<std::string> mass_matrix_path_arg(parser, "mass_matrix_path",
                                                       "File to save coarsened mass matrix to. If not set, the mass matrix is not computed", {'m',"mass_path"});
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (const args::ParseError& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  std::string filename = filename_arg ? args::get(filename_arg) : "../../meshes/spot.obj";
  int n_coarse_vertices = n_coarse_vertices_arg ? args::get(n_coarse_vertices_arg) : 500;
  double weight = weight_arg ? args::get(weight_arg) : 0;

  // load mesh
  MatrixXd VO;
  MatrixXi F, FO;
  {
    igl::read_triangle_mesh(filename, VO, FO);
  }

  if (n_coarse_vertices < 0 ) {
    std::cout << "Error: target number of vertices is negative: " << n_coarse_vertices << std::endl;
    exit(1);
  } else if (n_coarse_vertices >= VO.rows()) {
    std::cout << "Warning: target number of vertices is greater than input mesh size." << std::endl;
    std::cout <<"\t target number of vertices: " << n_coarse_vertices << std::endl;
    std::cout <<"\t input mesh size: " << VO.rows() << std::endl;
    n_coarse_vertices = VO.rows();
  }

  MatrixXi G;    // glue map
  MatrixXd l;    // edge lengths
  MatrixXd A;    // angular coordinates
  MatrixXi v2fs; // vertex to faceside map
  build_intrinsic_info(VO, FO, G, l, A, v2fs);
  F = FO;

  int total_removal = VO.rows() - n_coarse_vertices;
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
  for (Eigen::Index iV = 0; iV < V.rows(); iV++) {
    D(iV) = sin(V(iV, 0) * 10) * cos(V(iV, 1) * 5) + sin(V(iV, 2) * 8);
  }

  // prolongation
	VectorXd DO = P * D;

  if (prolongation_matrix_path_arg) {
    std::string path = args::get(prolongation_matrix_path_arg);
    saveMatrix(P, path);
  }

  if (laplace_matrix_path_arg) {
    std::string path = args::get(laplace_matrix_path_arg);
    SparseMatrix<double> L;
    igl::cotmatrix_intrinsic(l, F, L);
    saveMatrix(L, path);
  }

  if (mass_matrix_path_arg) {
    std::string path = args::get(mass_matrix_path_arg);
    SparseMatrix<double> M;
    igl::massmatrix_intrinsic(l, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    saveMatrix(M, path);
  }

  polyscope::init();
  polyscope::registerSurfaceMesh("input mesh", VO, FO);
  polyscope::getSurfaceMesh("input mesh")->addVertexScalarQuantity("prolonged rand function", DO);
  polyscope::registerSurfaceMesh("coarsened mesh (with wrong edge length)", V, F);
  polyscope::getSurfaceMesh("coarsened mesh (with wrong edge length)")->addVertexScalarQuantity("rand function", D);
  polyscope::view::lookAt(glm::vec3{1.5, 1.5, 3}, glm::vec3{0., 0., 0.});
  polyscope::show();
}

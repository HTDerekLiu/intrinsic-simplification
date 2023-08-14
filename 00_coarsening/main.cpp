#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "args/args.hxx"

#include <igl/read_triangle_mesh.h>
#include <igl/slice.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <build_intrinsic_info.h>
#include <chrono>
#include <coarsen_mesh.h>
#include <get_barycentric_points.h>
#include <remove_unreferenced_intrinsic.h>
#include <connected_components.h>

int main(int argc, char* argv[]) {
  using namespace Eigen;
  using namespace std;
  using namespace global_variables;
  using namespace std::chrono;

  // Configure the argument parser
  args::ArgumentParser parser("Intrinsic Coarsening");
  args::Positional<std::string> filename_arg(parser, "mesh",
                                              "Mesh to be coarsened. (default='../../meshes/bs_rest.obj')");
  args::Positional<int> n_coarse_vertices_arg(parser, "n_coarse_vertices",
                                        "number of coarse vertices to leave. (default='500')");
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

  std::string filename = filename_arg ? args::get(filename_arg) : "../../meshes/bs_rest.obj";
  int n_coarse_vertices = n_coarse_vertices_arg ? args::get(n_coarse_vertices_arg) : 500;

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

  // Check if mesh is connected
  VectorXi v_ids, f_ids;
  int n_components;
  connected_components(FO, G, n_components, v_ids, f_ids);
  if (n_components != 1) {
      std::cout << "WARNING: input mesh has " << n_components << " connected components. Simplification may behave unexpectedly when the input mesh is not connected." << std::endl;
  }

  std::cout << "starting simplification..." << std::endl;
  auto start = high_resolution_clock::now();

  int total_removal = VO.rows() - n_coarse_vertices;
  double weight = 0.0; // 0: pure curvature, 1: pure area
  MatrixXd BC;
  vector<vector<int>> F2V;
  coarsen_mesh(total_removal, weight, F, G, l, A, v2fs, BC, F2V);
  std::cout << "finished simplification..." << std::endl;

  auto stop     = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);

  cout << "removed " << total_removal << " vertices. runtime: " << ((float) duration.count() / 1000000.0) << " seconds" << endl;

  MatrixXd P;
  get_barycentric_points(VO, F, BC, F2V, P);

  // removed unreferenced vertices
  map<int, int> IMV, IMF;
  VectorXi vIdx, fIdx;
  remove_unreferenced_intrinsic(F,IMV, IMF, vIdx, fIdx);
  MatrixXd V;
  igl::slice(VO,vIdx,1,V);

  polyscope::init();
  auto psInputMesh = polyscope::registerSurfaceMesh("input mesh", VO, FO);
  if (n_components > 1) {
      VectorXd colors = VectorXd::Random(n_components);
      VectorXd vertex_colors;
      igl::slice(colors, v_ids, 1, vertex_colors);
      psInputMesh->addVertexScalarQuantity("component_color", vertex_colors);
  }
  polyscope::registerSurfaceMesh("coarsened mesh (with wrong edge length)", V, F);
  polyscope::registerPointCloud("barycentric points from the input", P);
  polyscope::view::lookAt(glm::vec3{1.5, 1.5, 3}, glm::vec3{0., 0., 0.});
  polyscope::show();
}

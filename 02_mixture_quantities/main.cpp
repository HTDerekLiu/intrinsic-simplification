#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "args/args.hxx"

#include <igl/read_triangle_mesh.h>

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
                                             "Mesh to be coarsened. (default='../../meshes/capsule.obj')");
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

  std::string filename = filename_arg ? args::get(filename_arg) : "../../meshes/capsule.obj";
  int n_coarse_vertices = n_coarse_vertices_arg ? args::get(n_coarse_vertices_arg) : 500;

  // load mesh
  MatrixXd VO;
  MatrixXi FO;
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

  int total_removal = VO.rows() - n_coarse_vertices;

  // pure curvature
  MatrixXd V_curvature;
  MatrixXi F_curvature = FO;
  {
    MatrixXi G;    // glue map
    MatrixXd l;    // edge lengths
    MatrixXd A;    // angular coordinates
    MatrixXi v2fs; // vertex to faceside map
    build_intrinsic_info(VO, FO, G, l, A, v2fs);

    // Check if mesh is connected
    VectorXi v_ids, f_ids;
    int n_components;
    connected_components(FO, G, n_components, v_ids, f_ids);
    if (n_components != 1) {
        std::cout << "WARNING: input mesh has " << n_components << " connected components. Simplification may behave unexpectedly when the input mesh is not connected." << std::endl;
    }

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
  // register results with polyscope
  auto psInput     = polyscope::registerSurfaceMesh("input mesh", VO, FO);
  auto psCurvature = polyscope::registerSurfaceMesh("curvature metric (with wrong edge length)", V_curvature, F_curvature);
  auto psMixture   = polyscope::registerSurfaceMesh("mixture metric (with wrong edge length)", V_mix, F_mix);
  auto psArea      =polyscope::registerSurfaceMesh("area metric (with wrong edge length)", V_area, F_area);
  polyscope::view::lookAt(glm::vec3{1.5, 1.5, 3}, glm::vec3{0., 0., 0.});

  // space out different results
  auto b_box = psInput->boundingBox();
  float width = std::get<1>(b_box).x - std::get<0>(b_box).x;
  // psInput->translate(     glm::vec3{-2.25, 0, 0} * width );
  // psCurvature->translate( glm::vec3{-.75 , 0, 0} * width );
  // psMixture->translate(   glm::vec3{ .75 , 0, 0} * width );
  // psArea->translate(      glm::vec3{ 2.25, 0, 0} * width );
  psInput->translate(     glm::vec3{-1.8, 0, 0} * width );
  psCurvature->translate( glm::vec3{-.6 , 0, 0} * width );
  psMixture->translate(   glm::vec3{ .6 , 0, 0} * width );
  psArea->translate(      glm::vec3{ 1.8, 0, 0} * width );

  polyscope::show();
}

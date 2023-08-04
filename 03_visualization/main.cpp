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
#include <remove_unreferenced_intrinsic.h>
#include <query_texture_barycentric.h>
#include <bake_texture.h>

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
  args::ValueFlag<int> texture_width_arg(parser, "texture_width",
                                     "texture width. (default='2048')", {'w',"texture_width"});
  args::ValueFlag<std::string> texture_path_arg(parser, "texture_path",
                                                "File to save texture to. Texture will be saved as a png. If not set, the texture is not saved", {'t',"texture_path"});
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
  int tex_width = texture_width_arg ? args::get(texture_width_arg) : 2048;

  // load mesh
  MatrixXd VO, UV, NV;
  MatrixXi F, FO, UF, NF;
  {
    igl::readOBJ(filename, VO, UV, NV, FO, UF, NF);
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

  int tex_size = tex_width*tex_width;
  MatrixXd bary_coords;
  VectorXi bary_faces;
  Matrix<bool, Dynamic, 1> hit_mask;
  query_texture_barycentric(UV, UF, tex_width, bary_faces, bary_coords, hit_mask);

  MatrixXi G;    // glue map
  MatrixXd l;    // edge lengths
  MatrixXd A;    // angular coordinates
  MatrixXi v2fs; // vertex to faceside map
  build_intrinsic_info(VO, FO, G, l, A, v2fs);
  F = FO;

  int total_removal = VO.rows() - n_coarse_vertices;
  MatrixXd BC;
  vector<vector<int>> F2V;

  int nV = v2fs.rows();
  int nF = F.rows();
  BC.resize(nV + tex_size, 3);
  BC.setConstant(0.0);
  F2V.resize(nF);

  for (int i = 0; i < tex_size; i++) {
    int idx = i + nV;
    if (hit_mask(i)) {
      BC.row(i + nV) = bary_coords.row(i);
      F2V[bary_faces(i)].push_back(idx);
    }
  }

  coarsen_mesh(total_removal, weight, F, G, l, A, v2fs, BC, F2V);

  cout << "removed " << total_removal << " vertices" << endl;

  std::vector<unsigned char> texture;
  bake_texture(texture, F, F2V, hit_mask, nV);

  if (texture_path_arg) {
    std::string texture_path = args::get(texture_path_arg);
    bake_texture(texture_path, texture);
  }

  // removed unreferenced vertices
  map<int, int> IMV, IMF;
  VectorXi vIdx, fIdx;
  remove_unreferenced_intrinsic(F, G, l, A, v2fs, F2V, IMV, IMF, vIdx, fIdx);
  MatrixXd V;
  igl::slice(VO,vIdx,1,V);

  // set up scene in polyscope
  polyscope::init();
  auto psMesh  = polyscope::registerSurfaceMesh("input mesh", VO, FO);
  polyscope::registerSurfaceMesh("coarsened mesh (with wrong edge length)", V, F);

  // convert parameterization to polyscope's desired input format
  Eigen::Matrix<glm::vec2, Dynamic, 1> parameterization(3 * UF.rows());
  for (int iF = 0; iF < UF.rows(); iF++) {
    for (int iC = 0; iC < 3; iC++) {
      parameterization(3 * iF + iC) = glm::vec2{UV(UF(iF, iC), 0), UV(UF(iF, iC), 1)};
    }
  }
  auto q = psMesh->addParameterizationQuantity("xz", parameterization)
                 ->setTexture(tex_width, tex_width, texture, polyscope::TextureFormat::RGBA8);
  q->setEnabled(true);
  q->setStyle(polyscope::ParamVizStyle::TEXTURE);
  q->setCheckerSize(1);

  polyscope::view::lookAt(glm::vec3{1.5, 1.5, 3}, glm::vec3{0., 0., 0.});
  polyscope::show();
}

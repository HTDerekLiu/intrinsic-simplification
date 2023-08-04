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
#include <mg_precompute.h>
#include <mass_matrix.h>
#include <cotan_Laplacian.h>
#include <min_quad_with_fixed_mg.h>

int main(int argc, char* argv[]) {
  using namespace Eigen;
  using namespace std;
  using namespace global_variables;
  using namespace std::chrono;

  // Configure the argument parser
  args::ArgumentParser parser("Intrinsic Multigrid Solver");
  args::Positional<std::string> filename_arg(parser, "mesh",
                                             "Mesh to be coarsened. (default='../../meshes/spot.obj')");
  args::ValueFlag<double> weight_arg(parser, "area_weight",
                                     "Influence of vertex area on coarsening. 0: none, 1: pure area weighting. (default='0.5')", {'w',"area_weight"});
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
  double weight = weight_arg ? args::get(weight_arg) : 0;

  // load mesh
  MatrixXd VO;
  MatrixXi F, FO;
  {
    igl::read_triangle_mesh(filename, VO, FO);
  }

  MatrixXi G;    // glue map
  MatrixXd l;    // edge lengths
  MatrixXd A;    // angular coordinates
  MatrixXi v2fs; // vertex to faceside map
  build_intrinsic_info(VO, FO, G, l, A, v2fs);
  F = FO;

  vector<mg_data> mg;
  double coarsening_ratio = 0.25;
  int approx_min_number_vertices = 500;
	mg_precompute(F,G,l,A,v2fs,coarsening_ratio,approx_min_number_vertices,weight,mg);

	// solve poisson problem
	SparseMatrix<double> L, M;
	cotan_Laplacian(F,l,L);
	mass_matrix(F,l,M);

	VectorXi b(3);
	b << 0, 1, 2;
	VectorXd bval(b.size());
	bval.setZero();

  // RHS
	VectorXd ones(VO.rows());
	ones.setOnes();
	VectorXd B = M * ones;
	for (int ii = 0; ii<b.size(); ii++)
		B(b(ii)) = bval(ii);

  // initial guess z0
	srand((unsigned int) time(0));
	VectorXd z0(VO.rows());
	z0.setRandom();

  // get some solver data
	VectorXd z = z0;
	min_quad_with_fixed_mg_data solverData;
	SimplicialLDLT<SparseMatrix<double>> coarseSolver;
	min_quad_with_fixed_mg_precompute(L,b,solverData, mg, coarseSolver);

	// multigrid solve
	vector<double> rHis;
	min_quad_with_fixed_mg_solve(solverData, B, bval, z0, coarseSolver, 1e-10, mg, z, rHis);

  polyscope::init();
	auto psMesh0 = polyscope::registerSurfaceMesh("mesh 0", VO, FO);
	psMesh0->addVertexScalarQuantity("solution", z);
	polyscope::show();
}

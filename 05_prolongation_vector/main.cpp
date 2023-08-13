#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "args/args.hxx"

#include <igl/read_triangle_mesh.h>
#include <igl/eigs.h>
#include <igl/massmatrix_intrinsic.h>
#include <igl/PI.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <build_intrinsic_info.h>
#include <chrono>
#include <coarsen_mesh.h>
#include <get_barycentric_points.h>
#include <remove_unreferenced_intrinsic.h>
#include <get_extrinsic_tangent_basis.h>
#include <get_vertex_connection_laplacian.h>
#include <get_vertex_vector_prolongation.h>
#include <save_complex_matrix.h>

// find eigenvector using simple power iteration
Eigen::VectorXcd complex_eigenvector(
    Eigen::SparseMatrix<std::complex<double>>& L,
    Eigen::SparseMatrix<std::complex<double>>& M)
{
    int N = L.rows();
    Eigen::VectorXcd u = Eigen::VectorXcd::Random(N);
    Eigen::VectorXcd x = u;

    auto residual = [&](const Eigen::VectorXcd& v) -> double {
        std::complex<double> lambda = v.dot(L * v);
        Eigen::VectorXcd err = L * v - lambda * M * v;
        return std::sqrt(std::abs(err.dot(M * err)));
    };

    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver(L);

    for (int i = 0; i < 500; i++) {
        // check for convergence
        if (residual(x) < 1e-9) break;

        // solve
        x = solver.solve(M * u);

        // normalize
        double norm = std::sqrt(std::abs(x.dot(M * x)));
        x /= norm;
        u = x;
    }

    return x;
}

int main(int argc, char* argv[]) {
    using namespace Eigen;
    using namespace std;
    using namespace global_variables;
    using namespace std::chrono;

    // Configure the argument parser
    args::ArgumentParser parser("Intrinsic Prolongation of Vector Fields");
    args::Positional<std::string> filename_arg(parser, "mesh",
                                              "Mesh to be coarsened. (default='../../meshes/dragon_fat.obj')");
    args::Positional<int> n_coarse_vertices_arg(parser, "n_coarse_vertices",
                                                "Number of coarse vertices to leave. (default='500')");
    args::ValueFlag<double> weight_arg(parser, "area_weight",
                                      "Influence of vertex area on coarsening. 0: none, 1: pure area weighting. (default='0')", {'w',"area_weight"});
    args::ValueFlag<std::string> prolongation_matrix_path_arg(parser, "prolongation_matrix_path",
                                      "File to save vector prolongation matrix to. If not set, the prolongation matrix is not saved", {'p',"prolongation_path"});
    args::ValueFlag<std::string> laplace_matrix_path_arg(parser, "laplace_matrix_path",
                                                              "File to save coarsened connection laplace matrix to. If not set, the laplace matrix is not saved", {'l',"laplace_path"});
    args::ValueFlag<std::string> mass_matrix_path_arg(parser, "mass_matrix_path",
                                                        "File to save coarsened vector mass matrix to. If not set, the mass matrix is not saved", {'m',"mass_path"});
    args::Flag no_viz_flag(parser, "no_viz", "Write requested output files without showing visualization", {'n', "no_viz"});
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

    std::string filename = filename_arg ? args::get(filename_arg) : "../../meshes/dragon_fat.obj";
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

    MatrixXi G, GO;        // glue map
    MatrixXd l, lO;        // edge lengths
    MatrixXd A, AO;        // angular coordinates
    MatrixXi v2fs, vO2fsO; // vertex to faceside map
    build_intrinsic_info(VO, FO, GO, lO, AO, vO2fsO);
    F = FO;
    G = GO;
    l = lO;
    A = AO;
    v2fs = vO2fsO;

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

    // get vector prolongation matrix
    SparseMatrix<std::complex<double>> P;
    get_vertex_vector_prolongation(FO, GO, lO, AO, vO2fsO,
                                   F,  G,  l,  A,  v2fs,
                                   vIdx, F2V, BC, P);

    // construct a hat function on the simplified mesh
    VectorXcd hat_function = VectorXd::Zero(V.rows());
    hat_function(0) = 1;

    // prolong the hat function to the fine mesh
    VectorXcd hatO = P * hat_function;

    // construct smooth vector field on simplified mesh using Globally-Optimal Direction Fields (Crane+ 2013)
    SparseMatrix<complex<double>> L;
    get_vertex_connection_laplacian(F, G, l, A, v2fs, L);

    SparseMatrix<double> Mreal;
    igl::massmatrix_intrinsic(l, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, Mreal);
    SparseMatrix<complex<double>> M = Mreal.cast<complex<double>>();

    VectorXcd f = complex_eigenvector(L, M);
    for (int iV = 0; iV < f.rows(); iV++) f(iV) /= std::abs(f(iV));

    // prolong the smooth vector field to the fine mesh
    VectorXcd fO = P * f;

    // get tangent bases to visualize fields
    Eigen::MatrixXd basisXO, basisYO, basisX, basisY;
    get_extrinsic_tangent_basis(VO, FO, AO, vO2fsO, basisXO, basisYO);
    get_extrinsic_tangent_basis(V, F, A, v2fs, basisX, basisY);

    if (prolongation_matrix_path_arg) {
        std::string path = args::get(prolongation_matrix_path_arg);
        save_complex_matrix(P, path);
    }

    if (laplace_matrix_path_arg) {
        std::string path = args::get(laplace_matrix_path_arg);
        save_complex_matrix(L, path);
    }

    if (mass_matrix_path_arg) {
        std::string path = args::get(mass_matrix_path_arg);
        save_complex_matrix(M, path);
    }

    if (no_viz_flag) {
        exit(0);
    }

    polyscope::init();
    auto psInput = polyscope::registerSurfaceMesh("input mesh", VO, FO);
    psInput->addVertexTangentVectorQuantity("prolonged hat function", hatO, basisXO, basisYO);
    psInput->addVertexTangentVectorQuantity("prolonged smooth field", fO, basisXO, basisYO)->setEnabled(true);
    auto psCoarse = polyscope::registerSurfaceMesh("coarsened mesh (with wrong edge length)", V, F);
    psCoarse->addVertexTangentVectorQuantity("hat function", hat_function, basisX, basisY);
    psCoarse->addVertexTangentVectorQuantity("smooth field", f, basisX, basisY)->setEnabled(true);
    psCoarse->setEnabled(false);
    polyscope::view::lookAt(glm::vec3{90, 90, 90}, glm::vec3{0., 30, 0.});
    polyscope::show();
}

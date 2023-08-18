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
#include <get_extrinsic_tangent_basis.h>
#include <get_vertex_connection_laplacian.h>
#include <get_prolongation.h>
#include <save_complex_matrix.h>
#include <get_vertex_vector_prolongation.h>
#include <save_matrix.h>
#include <connected_components.h>
#include <trace_geodesic.h>
#include <insert_degree_three_vertex.h>
#include <insert_ear_vertex.h>
#include <query_texture_barycentric.h>
#include <bake_texture.h>
#include <delaunay_refine.h>
#include <flip_to_delaunay.h>

enum class RefinementTime {Before = 0, After, Both, Neither};

int main(int argc, char* argv[]) {
    using namespace Eigen;
    using namespace std;
    using namespace global_variables;
    using namespace std::chrono;

    // Configure the argument parser
    args::ArgumentParser parser("Intrinsic simplification + Delaunay refinement");
    args::Positional<std::string> filename_arg(parser, "mesh", "Mesh to be coarsened."
                                               " (default='../../meshes/dragon_fat.obj')");
    args::Positional<int> n_coarse_vertices_arg(parser, "n_coarse_vertices",
                                                "Number of coarse vertices to leave."
                                                " (default='500')");
    args::ValueFlag<double> weight_arg(parser, "area_weight", "Influence of vertex area"
                                       " on coarsening. 0: none, 1: pure area weighting."
                                       " (default='0')", {'w',"area_weight"});
    args::ValueFlag<std::string> refinement_time_arg(parser, "refinement_time", "Whether to"
                                                     " perform intrinsic Delaunay refinement"
                                                     " before simplification or after"
                                                     " simplification. (Options: BEFORE, AFTER,"
                                                     " BOTH, NEITHER; default=AFTER)",
                                                     {'r', "refinement_time"}
                                                     );
    args::ValueFlag<int> texture_width_arg(parser, "texture_width", "Texture width. Set to"
                                           " -1 to disable texture visualization"
                                           " (default='2048')", {'u',"texture_width"});
    args::ValueFlag<std::string> texture_path_arg(parser, "texture_path", "File to save texture"
                                                  " to. Texture will be saved as a png. If not"
                                                  " set, the texture is not saved",
                                                  {'t',"texture_path"});
    args::ValueFlag<std::string> prolongation_matrix_path_arg(parser, "prolongation_matrix_path",
                                      "File to save prolongation matrix""to. If not set, the"
                                      " prolongation matrix is not saved",
                                      {'p',"prolongation_path"});
    args::ValueFlag<std::string> laplace_matrix_path_arg(parser, "laplace_matrix_path",
                                      "File to save coarsened laplace matrix to. If not set, the"
                                      " laplace matrix is not computed", {'l',"laplace_path"});
    args::ValueFlag<std::string> mass_matrix_path_arg(parser, "mass_matrix_path",
                                      "File to save coarsened mass matrix to. If not set, the"
                                      " mass matrix is not computed", {'m',"mass_path"});
    args::ValueFlag<std::string> v_prolongation_matrix_path_arg(parser,
                                      "vector_prolongation_matrix_path",
                                      "File to save vector prolongation matrix to. If not set,"
                                      " the prolongation matrix is not saved",
                                      {"vp","vector_prolongation_path"});
    args::ValueFlag<std::string> v_laplace_matrix_path_arg(parser,
                                      "connection_laplace_matrix_path",
                                      "File to save coarsened connection laplace matrix to."
                                      " If not set, the laplace matrix is not computed",
                                      {"cl","connection_laplace_path"});
    args::ValueFlag<std::string> v_mass_matrix_path_arg(parser, "vector_mass_matrix_path",
                                     "File to save coarsened vector mass matrix to. If not set,"
                                     " the mass matrix is not computed",
                                     {"vm", "vector_mass_path"});
    args::Flag no_viz_flag(parser, "no_viz", "Write requested output files without showing"
                           " visualization", {'n', "no_viz"});
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

    std::string filename = filename_arg ? args::get(filename_arg)
                                        : "../../meshes/dragon_fat.obj";
    int n_coarse_vertices = n_coarse_vertices_arg ? args::get(n_coarse_vertices_arg) : 500;
    double weight = weight_arg ? args::get(weight_arg) : 0;
    int tex_width = texture_width_arg ? args::get(texture_width_arg) : 2048;

    RefinementTime refinement_time = RefinementTime::After;
    if (refinement_time_arg) {
        std::string refinement_time_str = args::get(refinement_time_arg);
        if (refinement_time_str == "BEFORE") {
            refinement_time = RefinementTime::Before;
        } else if (refinement_time_str == "AFTER") {
            refinement_time = RefinementTime::After;
        } else if (refinement_time_str == "BOTH") {
            refinement_time = RefinementTime::Both;
        } else if (refinement_time_str == "NEITHER") {
            refinement_time = RefinementTime::Neither;
        } else {
            std::cout << "Error: unrecognized refinement time '" << refinement_time_str
                      << "'" << std::endl;
            exit(1);
        }
    }

    // load mesh
    MatrixXd VO, UV, NV;
    MatrixXi F, FO, UF, NF;
    {
        igl::readOBJ(filename, VO, UV, NV, FO, UF, NF);

    }
    if (UV.rows() == 0) {
        std::cout << "Error: input mesh has no UV coordinates" << std::endl;
        exit(1);
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

    bool using_texture = tex_width > 0;

    // ======================================================
    //                    initialization
    // ======================================================
    int tex_size = (using_texture) ? tex_width * tex_width : 0;
    MatrixXd bary_coords;
    VectorXi bary_faces;
    Matrix<bool, Dynamic, 1> hit_mask;
    if (using_texture) {
        query_texture_barycentric(UV, UF, tex_width, bary_faces, bary_coords, hit_mask);
    } else {
        std::cout << "Skipping texture generation since width was set to "
                  << tex_width << std::endl;
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

    // Check if mesh is connected
    VectorXi v_ids, f_ids;
    int n_components;
    connected_components(FO, G, n_components, v_ids, f_ids);
    if (n_components != 1) {
        std::cout << "WARNING: input mesh has " << n_components << " connected components."
          " Simplification may behave unexpectedly when the input mesh is not connected."
                  << std::endl;
    }

    MatrixXd BC;
    vector<vector<int>> F2V;

    int nVO = v2fs.rows();
    int nV = nVO;
    int nF = F.rows();
    BC.resize(nVO + tex_size, 3);
    BC.setConstant(DOUBLE_INF);
    F2V.resize(nF);

    size_t iP = 0;
    for (int i = 0; i < tex_size; i++) {
        if (hit_mask(i)) {
            BC.row(nVO + iP) = bary_coords.row(i);
            F2V[bary_faces(i)].push_back(nVO + iP);
            iP++;
        }
    }
    BC.conservativeResize(nVO + iP, 3); // only keep active pixels

    // ======================================================
    //                    retriangulation
    // ======================================================

    // == refinement before coarsening
    if (refinement_time == RefinementTime::Before
        || refinement_time == RefinementTime::Both) {
        cout << "refining before coarsening ... " << endl;
        int n_insertions = delaunay_refine(F, G, l, A, v2fs, BC, F2V);
        cout << "   inserted " << n_insertions << " new vertices" << endl;

        nV = 0; // recount number of vertices after refinement
        for (int iV = 0; iV < v2fs.rows(); iV++) {
            // filter out deleted vertices
            if (v2fs(iV, 0) != GHOST_INDEX) nV++;
        }

        // rebuild BC array to make space for new vertices
        MatrixXd new_BC;
        int index_shift = v2fs.rows() - nVO;
        new_BC.resize(v2fs.rows() + iP, 3);
        new_BC.setConstant(DOUBLE_INF);

        for (int iF = 0; iF < F.rows(); iF++) {
          // skip deleted faces
          if (F(iF, 0) == GHOST_INDEX) continue;

          for (int& iP : F2V[iF]) {
            new_BC.row(iP + index_shift) = BC.row(iP);
            iP += index_shift;
          }
        }

        BC = new_BC;
        nV = v2fs.rows();
    }

    // == coarsening
    int total_removal = nV - n_coarse_vertices;
    cout << "coarsening ... " << endl;
    coarsen_mesh(total_removal, weight, F, G, l, A, v2fs, BC, F2V);

    // == refinement after coarsening
    if (refinement_time == RefinementTime::After
        || refinement_time == RefinementTime::Both) {
        cout << "refining after coarsening ... " << endl;
        int n_insertions = delaunay_refine(F, G, l, A, v2fs, BC, F2V);
        cout << "   inserted " << n_insertions << " new vertices" << endl;
    }

    // ======================================================
    //               compute assorted outputs
    // ======================================================
    std::vector<unsigned char> texture;
    if (using_texture) {
        bake_texture(texture, F, F2V, hit_mask, nV);

        if (texture_path_arg) {
            std::string texture_path = args::get(texture_path_arg);
            bake_texture(texture_path, texture);
        }
    }

    // removed unreferenced vertices
    map<int, int> IMV, IMF;
    VectorXi vIdx, fIdx;
    remove_unreferenced_intrinsic(F, G, l, A, v2fs, F2V, IMV, IMF, vIdx, fIdx);

    // get scalar prolongation matrix
    SparseMatrix<double> P;
    get_prolongation(F, BC.block(0, 0, nVO, 3), F2V, IMV, P);

    // create a scalar hat function on the simplified mesh
    // int i_hat = (int)(vIdx.rows() * 0.45);
    int i_hat = (int)(vIdx.rows() * 0.25);
    VectorXd hat_function = VectorXd::Zero(vIdx.rows());
    hat_function(i_hat) = 1;

    // prolong the hat function to the fine mesh
    VectorXd hatO = P * hat_function;

    // get vector prolongation matrix
    SparseMatrix<std::complex<double>> vP;
    get_vertex_vector_prolongation(FO, GO, lO, AO, vO2fsO,
                                   F,  G,  l,  A,  v2fs,
                                   vIdx, F2V, BC.block(0, 0, nVO, 3), vP);

    // construct a vector hat function on the simplified mesh
    VectorXcd vector_hat_function = VectorXd::Zero(vIdx.rows());
    vector_hat_function(i_hat) = 1;

    // prolong the hat function to the fine mesh
    VectorXcd vector_hatO = vP * vector_hat_function;

    // get a tangent basis to visualize fields
    Eigen::MatrixXd basisXO, basisYO;
    get_extrinsic_tangent_basis(VO, FO, AO, vO2fsO, basisXO, basisYO);

    if (prolongation_matrix_path_arg) {
        std::string path = args::get(prolongation_matrix_path_arg);
        save_matrix(P, path);
    }

    if (laplace_matrix_path_arg) {
        std::string path = args::get(laplace_matrix_path_arg);
        SparseMatrix<double> L;
        igl::cotmatrix_intrinsic(l, F, L);
        save_matrix(L, path);
    }

    if (mass_matrix_path_arg) {
        std::string path = args::get(mass_matrix_path_arg);
        SparseMatrix<double> M;
        igl::massmatrix_intrinsic(l, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
        save_matrix(M, path);
    }

    if (v_prolongation_matrix_path_arg) {
        std::string path = args::get(prolongation_matrix_path_arg);
        save_complex_matrix(vP, path);
    }

    if (v_laplace_matrix_path_arg) {
      std::string path = args::get(v_laplace_matrix_path_arg);
      SparseMatrix<complex<double>> cL;
      get_vertex_connection_laplacian(F, G, l, A, v2fs, cL);
      save_complex_matrix(cL, path);
    }

    if (v_mass_matrix_path_arg) {
      std::string path = args::get(v_mass_matrix_path_arg);
      SparseMatrix<double> Mreal;
      igl::massmatrix_intrinsic(l, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, Mreal);
      SparseMatrix<complex<double>> vM = Mreal.cast<complex<double>>();
      save_complex_matrix(vM, path);
    }

    if (no_viz_flag) {
        exit(0);
    }

    polyscope::init();
    auto psMesh  = polyscope::registerSurfaceMesh("input mesh", VO, FO);

    if (using_texture) {
        // convert parameterization to polyscope's desired input format
        Eigen::Matrix<glm::vec2, Dynamic, 1> parameterization(3 * UF.rows());
        for (int iF = 0; iF < UF.rows(); iF++) {
            for (int iC = 0; iC < 3; iC++) {
                parameterization(3 * iF + iC) = glm::vec2{UV(UF(iF, iC), 0), UV(UF(iF, iC), 1)};
            }
        }
        auto q = psMesh->addParameterizationQuantity("intrinsic triangulation",
                                                     parameterization)
                  ->setTexture(tex_width, tex_width, texture, polyscope::TextureFormat::RGBA8);
        q->setEnabled(true);
        q->setStyle(polyscope::ParamVizStyle::TEXTURE);
        q->setCheckerSize(1);
    }

    psMesh->addVertexScalarQuantity("prolonged scalar hat function", hatO);
    psMesh->addVertexTangentVectorQuantity("prolonged vector hat function",
                                           hatO, basisXO, basisYO);
    polyscope::show();
}

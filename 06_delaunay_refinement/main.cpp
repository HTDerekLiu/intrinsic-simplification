#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "args/args.hxx"

#include <igl/read_triangle_mesh.h>
#include <igl/eigs.h>
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
                                                     " BOTH, NEITHER; default=BEFORE)",
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
                                      "File to save vector prolongation matrix to. If not set,"
                                      " the prolongation matrix is not saved",
                                      {'p',"prolongation_path"});
    args::ValueFlag<std::string> laplace_matrix_path_arg(parser, "laplace_matrix_path",
                                      "File to save coarsened connection laplace matrix to."
                                      " If not set, the laplace matrix is not saved",
                                      {'l',"laplace_path"});
    args::ValueFlag<std::string> mass_matrix_path_arg(parser, "mass_matrix_path",
                                     "File to save coarsened vector mass matrix to. If not set,"
                                     " the mass matrix is not saved", {'m',"mass_path"});
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

    std::string filename = filename_arg ? args::get(filename_arg) : "../../meshes/spot.obj";
    int n_coarse_vertices = n_coarse_vertices_arg ? args::get(n_coarse_vertices_arg) : 500;
    double weight = weight_arg ? args::get(weight_arg) : 0;
    int tex_width = texture_width_arg ? args::get(texture_width_arg) : 2048;

    RefinementTime refinement_time = RefinementTime::Before;
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

    // check if face has been deleted
    auto is_dead_face = [&](int iF) -> bool { return F(iF, 0) == GHOST_INDEX; };

    // Check if mesh is connected
    VectorXi v_ids, f_ids;
    int n_components;
    connected_components(FO, G, n_components, v_ids, f_ids);
    if (n_components != 1) {
        std::cout << "WARNING: input mesh has " << n_components << " connected components."
          " Simplification may behave unexpectedly when the input mesh is not connected."
                  << std::endl;
    }

    // test insert_degree_three_vertex
    if (false) {
        auto bc_to_position = [](int f,
                                 const Eigen::Vector3d & bc,
                                 const Eigen::MatrixXd & V,
                                 const Eigen::MatrixXi & F) -> Eigen::Vector3d {
            return bc(0) * V.row(F(f, 0)) + bc(1) * V.row(F(f, 1)) + bc(2) * V.row(F(f, 2));
        };

        auto inside_triangle = [](const Vector3d& v) -> bool {
          return v(0) >= 0 && v(1) >= 0 && v(2) >= 0;
        };

        polyscope::init();
        auto psInput = polyscope::registerSurfaceMesh("input mesh", VO, FO);
        int f = 3609; // spot.obj
        std::cout << "Face vertices: " << F.row(f) << std::endl;
        Eigen::Vector3d b{0.1, 0.3, 0.6};
        Eigen::Vector3d x = bc_to_position(f, b, VO, FO);
        auto psNewPoint = polyscope::registerPointCloud("inserted point", x.transpose());

        int nPts = 25;
        Eigen::MatrixXd BC(nPts, 3), pt_positions(nPts, 3);
        std::vector<std::vector<int>> F2V(F.rows(), std::vector<int>{});
        for (int iP = 0; iP < nPts; iP++) {
            Eigen::Vector3d pb = Eigen::Vector3d::Random().cwiseAbs();
            pb /= pb.sum();
            BC.row(iP) = pb;
            F2V[f].push_back(iP);
            pt_positions.row(iP) = bc_to_position(f, pb, VO, FO);
        }
        double r = 0.001;
        polyscope::registerPointCloud("original points", pt_positions)->setPointRadius(r);

        int v = insert_degree_three_vertex(f, b, F, G, l, A, v2fs, BC, F2V);
        Eigen::MatrixXd V = VO;
        V.conservativeResize(V.rows()+1, 3);
        V.row(v) = x;

        auto psSubdivided = polyscope::registerSurfaceMesh("subdivided mesh", V, F);
        Eigen::MatrixXd aa(F2V[f].size(), 3), bb(F2V[F.rows()-2].size(), 3),
                        cc(F2V[F.rows()-1].size(), 3);
        int ii = 0;
        for (int iP : F2V[f]) {
          aa.row(ii++) = bc_to_position(f, BC.row(iP), V, F);
        }
        ii = 0;
        for (int iP : F2V[F.rows()-2]) {
          bb.row(ii++) = bc_to_position(F.rows()-2, BC.row(iP), V, F);
        }
        ii = 0;
        for (int iP : F2V[F.rows()-1]) {
          cc.row(ii++) = bc_to_position(F.rows()-1, BC.row(iP), V, F);
        }
        polyscope::registerPointCloud("subdivided points a", aa)->setPointRadius(r);
        polyscope::registerPointCloud("subdivided points b", bb)->setPointRadius(r);
        polyscope::registerPointCloud("subdivided points c", cc)->setPointRadius(r);

        polyscope::show();
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
    int nPix = BC.rows();

    // pre-refinement
    if (refinement_time == RefinementTime::Before
        || refinement_time == RefinementTime::Both) {
        delaunay_refine(F, G, l, A, v2fs, BC, F2V);

        nV = 0; // recount number of vertices after refinement
        for (int iV = 0; iV < v2fs.rows(); iV++) {
            // filter out deleted vertices
            if (v2fs(iV, 0) != GHOST_INDEX) nV++;
        }


        for (int f = 0; f < F.rows(); f++) {
            if (is_dead_face(f) && !F2V[f].empty())
              throw std::runtime_error("dead face contains pts");
        }

        // rebuild BC array to make space for new vertices
        MatrixXd new_BC;
        int index_shift = v2fs.rows() - nVO;
        new_BC.resize(v2fs.rows() + iP, 3);
        new_BC.setConstant(DOUBLE_INF);

        for (int iF = 0; iF < F.rows(); iF++) {
          if (is_dead_face(iF)) continue;

          for (int& iP : F2V[iF]) {
            new_BC.row(iP + index_shift) = BC.row(iP);
            iP += index_shift;
          }
        }

        BC = new_BC;
        nVO = v2fs.rows();
    }

    polyscope::init();
    auto psMesh  = polyscope::registerSurfaceMesh("input mesh", VO, FO);

    auto snapshot = [&](std::string name, int padded_length = 3) {
      cout << "baking texture with nVO = " << nVO
           << ", iP = " << iP
           << ", F2V.size() = " << F2V.size()
           << ", texture size = " << hit_mask.size() << endl;
      std::vector<unsigned char> texture;
      bake_texture(texture, F, F2V, hit_mask, nVO);

      // convert parameterization to polyscope's desired input format
      Eigen::Matrix<glm::vec2, Dynamic, 1> parameterization(3 * UF.rows());
      for (int iF = 0; iF < UF.rows(); iF++) {
        for (int iC = 0; iC < 3; iC++) {
          parameterization(3 * iF + iC) = glm::vec2{UV(UF(iF, iC), 0), UV(UF(iF, iC), 1)};
        }
      }
      name.insert(name.begin(), padded_length - name.size(), ' ');
      auto q = psMesh->addParameterizationQuantity("intrinsic triangulation " + name, parameterization)
        ->setTexture(tex_width, tex_width, texture, polyscope::TextureFormat::RGBA8);
      q->setEnabled(true);
      q->setStyle(polyscope::ParamVizStyle::TEXTURE);
      q->setCheckerSize(1);

      // polyscope::registerSurfaceMesh2D("UV " + name, UV, UF);
    };

    if (true) {
        // int step = 1;
        // int f_watch = 1371;
        // for (int i = 0; i < 5; i+=step) {
        //   if (!is_dead_face(f_watch)){
        //     cout << "iteration " << i << " | face " << f_watch << " contains points: " << endl;
        //     if (f_watch >= F2V.size()) {
        //       std::cout << "f_watch " << f_watch << " >= F2V.size() " << F2V.size() << endl;
        //     }
        //     for (int iP : F2V[f_watch]) {
        //       cout << "  " << iP <<  " (" << (iP - nVO) << ") | " << BC.row(iP) << endl;
        //     }}
        //     snapshot(std::to_string(i));

        //     // check that we've kept all points
        //     int nPts = 0;
        //     for (int iF = 0; iF < F.rows(); iF++) {
        //         if (F(iF, 0) < 0) continue;
        //         for (int idx : F2V[iF]) {
        //           if (idx >= nVO) {
        //             nPts++;
        //           }
        //         }
        //     }
        //     if (nPts + nVO != BC.rows()) {
        //       std::cout << "Error, found " << nPts << " points, but there should be " << (BC.rows() - nVO) << std::endl;
        //     }
        //     cout << "started with " << nPix << " bary coords, and now we have " << BC.rows() << endl;
        //     cout << "done with info dump" << endl << endl;

        //     coarsen_mesh(step, weight, F, G, l, A, v2fs, BC, F2V);
        // }
        // coarsening
        int total_removal = nV - n_coarse_vertices;
        coarsen_mesh(total_removal, weight, F, G, l, A, v2fs, BC, F2V);

        // cout << "removed " << total_removal << " vertices" << endl;
    }

    // post-refinement
    if (refinement_time == RefinementTime::After
        || refinement_time == RefinementTime::Both) {
        delaunay_refine(F, G, l, A, v2fs, BC, F2V);
    }

    // removed unreferenced vertices
    // map<int, int> IMV, IMF;
    // VectorXi vIdx, fIdx;
    // remove_unreferenced_intrinsic(F, G, l, A, v2fs, F2V, IMV, IMF, vIdx, fIdx);
    // MatrixXd V; // TODO: interpolate vertex positions to inserted vertices
    // igl::slice(VO,vIdx,1,V);


    if (using_texture) {
        std::vector<unsigned char> texture;
        bake_texture(texture, F, F2V, hit_mask, nVO);

        if (texture_path_arg) {
            std::string texture_path = args::get(texture_path_arg);
            bake_texture(texture_path, texture);
        }

        // convert parameterization to polyscope's desired input format
        Eigen::Matrix<glm::vec2, Dynamic, 1> parameterization(3 * UF.rows());
        for (int iF = 0; iF < UF.rows(); iF++) {
            for (int iC = 0; iC < 3; iC++) {
                parameterization(3 * iF + iC) = glm::vec2{UV(UF(iF, iC), 0), UV(UF(iF, iC), 1)};
            }
        }
        auto q = psMesh->addParameterizationQuantity("intrinsic triangulation", parameterization)
                  ->setTexture(tex_width, tex_width, texture, polyscope::TextureFormat::RGBA8);
        q->setEnabled(true);
        q->setStyle(polyscope::ParamVizStyle::TEXTURE);
        q->setCheckerSize(1);

        // polyscope::registerSurfaceMesh2D("uv", UV, UF);
    }

    polyscope::show();

    /*
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

    // // construct smooth vector field on simplified mesh using Globally-Optimal Direction Fields (Crane+ 2013)
    SparseMatrix<complex<double>> L;
    get_vertex_connection_laplacian(F, G, l, A, v2fs, L);

    SparseMatrix<double> Mreal;
    igl::massmatrix_intrinsic(l, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, Mreal);
    SparseMatrix<complex<double>> M = Mreal.cast<complex<double>>();

    // VectorXcd f = complex_eigenvector(L, M);
    // for (int iV = 0; iV < f.rows(); iV++) f(iV) /= std::abs(f(iV));

    // // prolong the smooth vector field to the fine mesh
    // VectorXcd fO = P * f;

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
    // psInput->addVertexTangentVectorQuantity("prolonged smooth field", fO, basisXO, basisYO)->setEnabled(true);
    auto psCoarse = polyscope::registerSurfaceMesh("coarsened mesh (with wrong edge length)", V, F);
    psCoarse->addVertexTangentVectorQuantity("hat function", hat_function, basisX, basisY);
    // psCoarse->addVertexTangentVectorQuantity("smooth field", f, basisX, basisY)->setEnabled(true);
    psCoarse->setEnabled(false);
    polyscope::view::lookAt(glm::vec3{90, 90, 90}, glm::vec3{0., 30, 0.});
    polyscope::show();
 */
}

#include "query_texture_barycentric.h"
void query_texture_barycentric(
    const Eigen::MatrixXd & UV,
    const Eigen::MatrixXi & UF,
    const int & tex_width,
    Eigen::VectorXi & bary_faces,
    Eigen::MatrixXd & bary_coords,
    Eigen::Matrix<bool, Eigen::Dynamic, 1> & hit_mask)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    igl::AABB<MatrixXd, 2> tree;
    tree.init(UV, UF);

    int tex_size = tex_width * tex_width;

    bary_coords.resize(tex_size, 3);
    bary_faces.resize(tex_size);
    bary_faces.setConstant(GHOST_INDEX);
    hit_mask.resize(tex_size);

    for (int j = 0; j < tex_width; j++) {
        for (int i = 0; i < tex_width; i++) {
            int coord = j*tex_width + i;

            double u = ((double)i + 0.5) / (double)tex_width;
            double v = ((double)j + 0.5) / (double)tex_width;
            RowVector2d uv = RowVector2d(u, v);

            // find which face it is on
            hit_mask(coord) = false;
            int f;
            RowVector2d C;
            double dst = tree.squared_distance(UV, UF, uv, f, C);
            if (dst < 1e-4) {
                bool is_degenerated;
                bool is_inside;
                double distance_to_valid;
                Vector3d bc = compute_barycentric_robust(C, UV.row(UF(f, 0)), UV.row(UF(f, 1)), UV.row(UF(f, 2)), 1e-7, is_degenerated, is_inside,distance_to_valid);
                if (!is_degenerated && (distance_to_valid<1e-4)) // make sure to only append valid barycentric
                {
                    bary_coords.row(coord) = bc;
                    bary_faces(coord) = f;
                    hit_mask(coord) = true;
                }
            }
        }
    }
}

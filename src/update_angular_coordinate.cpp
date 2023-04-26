#include "update_angular_coordinate.h"
void update_angular_coordinate(
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::MatrixXd & A)
{
    using namespace std;
    using namespace Eigen;

    bool is_boundary_v;
    vector<Vector2i> fs_list;
    vertex_one_ring_face_sides(G, v2fs, v, fs_list, is_boundary_v);

    // compute real angles
    VectorXd angles(fs_list.size());
    for (int ii=0; ii<fs_list.size(); ii++)
    {
        angles(ii) = opposite_corner_angle(l, next(fs_list[ii]));
        if (isnan(angles(ii))) // if invalid triangles, change notion
            return;
    }

    // normalize angles
    if (is_boundary_v)
        angles = angles / (angles.sum() + M_PI) * 2 * M_PI;
    else
        angles = angles / angles.sum() * 2 * M_PI;

    
    // put angular coorinates back to A
    int num_corners = angles.size();
    VectorXd angles_cumsum(num_corners);
    igl::cumsum(angles, 1, angles_cumsum);

    // keep A_vi_cumsum
    VectorXd angles_coord(num_corners);
    angles_coord << 0.0, angles_cumsum.head(num_corners-1);

    // get Amin
    Vector2i fs_min = v2fs.row(v).transpose();
    double angle_min = A(fs_min(0), fs_min(1));

    // put it back to A
    for (auto it = std::begin(fs_list); it != std::end(fs_list); ++it)
    {
        Vector2i fs = *it;
        int idx = it - fs_list.begin();
        A(fs(0), fs(1)) = angles_coord(idx) + angle_min;
    }

}
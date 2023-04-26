#include "build_angular_coordinates.h"

void build_angular_coordinates(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs)
{
    using namespace Eigen;
    using namespace std;

    int nF = F.rows();
    int nV = F.maxCoeff() + 1;

    // resize A and v2fs
    A.resize(nF,3);
    v2fs.resize(nV,2);

    // an array to store which vertices have been visited
    Array<bool,Dynamic,1> visited;
    visited.setConstant(nV,false); 

    for (int f=0; f<nF; f++)
    {
        for (int s=0; s<3; s++)
        {
            // this face side
            Vector2i fs;
            fs << f, s;

            // this vertex
            int vi = F(f,s);

            // if not visited
            if (visited(vi) == false)
            {
                visited(vi) = true; // set vi to visited

                // get list of one-ring face sides
                std::vector<Vector2i> fs_list;
                bool is_boundary_vertex;
                {
                    vertex_one_ring_face_sides(G,fs,fs_list, is_boundary_vertex);
                }

                // compute angles for each corner attached to vi
                int num_corners = fs_list.size();
                VectorXd A_vi(num_corners);
                {
                    for (auto it = std::begin(fs_list); it != std::end(fs_list); ++it)
                    {
                        Vector2i fs = *it;
                        int idx = it - fs_list.begin();

                        // compute opposite corner angle
                        fs = next(fs);
                        double angle = opposite_corner_angle(l, fs);

                        // assign back to Avi(idx)
                        A_vi(idx) = angle;
                    }

                    // normalize the angular coordinates
                    if (is_boundary_vertex) 
                        A_vi = A_vi / (A_vi.sum()+M_PI) * M_PI*2;
                    else
                        A_vi = A_vi / (A_vi.sum()) * M_PI*2;
                }

                // put angular coorinates back to A
                {
                    VectorXd A_vi_cumsum(num_corners);
                    igl::cumsum(A_vi, 1, A_vi_cumsum);

                    // keep A_vi_cumsum
                    VectorXd A_coord(num_corners);
                    A_coord << 0.0, A_vi_cumsum.head(num_corners-1);

                    // put it back to A
                    for (auto it = std::begin(fs_list); it != std::end(fs_list); ++it)
                    {
                        Vector2i fs = *it;
                        int idx = it - fs_list.begin();

                        A(fs(0), fs(1)) = A_coord(idx);
                    }
                }

                // assign v2fs
                Vector2i fs = fs_list[0];
                v2fs.row(vi) = fs.transpose();
            }
        }
    }
}
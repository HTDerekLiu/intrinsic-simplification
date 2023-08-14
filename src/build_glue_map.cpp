#include "build_glue_map.h"

void build_glue_map(
    const Eigen::MatrixXi & F,
    Eigen::MatrixXi & G)
{
    using namespace Eigen;
    using namespace global_variables;

    int nF = F.rows(); // number of faces

    // Build a temparary list S of all face-sides, each row of S contains(vi,vj,f,s)
    MatrixXi S;
    {
        int nS = 3 * nF; // number of sides

        MatrixXi S_unsorted;
        S_unsorted.resize(nS, 4);
        S_unsorted.fill(0);

        for (int f=0; f<nF; f++)
        {
            for (int s=0; s<3; s++)
            {
                // initialize face side
                Vector2i fs;
                fs << f, s;

                // get vertices of this face side
                int vi, vj;
                get_face_side_vertices(F,fs,vi,vj);

                // assign to S
                S_unsorted(f*3+s, 0) = std::min(vi, vj);
                S_unsorted(f*3+s, 1) = std::max(vi, vj);
                S_unsorted(f*3+s, 2) = f;
                S_unsorted(f*3+s, 3) = s;
            }
        }
        igl::sortrows(S_unsorted, true, S);
    }

      // Build the gluing map G
    G.resize(nF, 3*2);
    G.fill(0);
    {
        int ii = 0;
        int nS = S.rows();
        Vector2i fs0, fs1;
        while (ii < nS)
        {
            if (ii == (nS-1))
            {
                // last row (must be a boundary face side)
                fs0 << S(ii  , 2), S(ii  , 3);
                glue_face_sides(fs0, GHOST_FACE_SIDE, G);
                ii += 1;
            }
            else if (S(ii, 0) != S(ii+1, 0) || S(ii, 1) != S(ii+1, 1))
            {
                // boundary face side
                fs0 << S(ii  , 2), S(ii  , 3);
                glue_face_sides(fs0, GHOST_FACE_SIDE, G);
                ii += 1;
            }
            else if (ii + 2 < nS && S(ii, 0) == S(ii+2, 0) && S(ii, 1) == S(ii+2, 1))
            {
                // nonmanifold face side
                throw std::invalid_argument("Input mesh is not manifold. Three faces were found"
                                            "along the edge between vertices "
                                            + std::to_string(S(ii, 0))
                                            + " and " + std::to_string(S(ii, 1)));
            }
            else
            {
                // interior face side
                fs0 << S(ii  , 2), S(ii  , 3);
                fs1 << S(ii+1, 2), S(ii+1, 3);

                // check orientation
                int vi, vj, wi, wj;
                get_face_side_vertices(F,fs0,vi,vj);
                get_face_side_vertices(F,fs1,wi,wj);
                if (vi == wi && vj == wj) {
                    throw std::invalid_argument("Input mesh is not oriented properly."
                                                "Two edges were found going from vertex "
                                                + std::to_string(vi) + " to vertex "
                                                + std::to_string(vj));
                }

                glue_face_sides(fs0, fs1, G);
                ii += 2;
            }
        }
    }
}

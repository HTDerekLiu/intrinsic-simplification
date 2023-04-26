#include "face_n_ring_faces.h"
void face_n_ring_faces(
    const Eigen::MatrixXi & G,
    const int & f_src,
    const int & num_rings,
    Eigen::VectorXi & fIdx)
{
    using namespace std;
    using namespace Eigen; 
    using namespace global_variables;

    // create local face sides
    set<int> fIdx_set;
    fIdx_set.insert(f_src);

    for (int iter=0; iter<num_rings; iter++)
    {
        set<int> to_append_fIdx_set; 
        for (int f : fIdx_set)
        {
            for (int s=0; s<3; s++)
            {
                Vector2i fs, fs_twin;
                fs << f, s;
                fs_twin = twin(G, fs);
                if (!is_same_face_side(fs_twin, GHOST_FACE_SIDE))
                    to_append_fIdx_set.insert(fs_twin(0)); // insert to a set (will avoid duplication)
            }
        }
        fIdx_set.insert(to_append_fIdx_set.begin(), to_append_fIdx_set.end());
    }

    // turn set into vectorXi
    fIdx.resize(fIdx_set.size());
    int ii = 0;
    for (int f : fIdx_set)
    {
        fIdx(ii) = f;
        ii += 1;
    }
}
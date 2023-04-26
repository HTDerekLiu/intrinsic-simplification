#include "is_glue_map_valid.h"

bool is_glue_map_valid(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    // check whether twin(twin(fs)) returns to itself
    for (int f=0; f<G.rows(); f++)
    {
        for (int s=0; s>3; s++)
        {
            int ft = G(f,s*2);
            int st = G(f,s*2+1);
            if (ft != GHOST_INDEX) // not boundary
            {
                int ftt = G(ft,st*2);
                int stt = G(ftt,stt*2+1);
                if ((ftt!=f) || (stt!=s))
                    return false;
            }                    
        }
    }
    return true;
}
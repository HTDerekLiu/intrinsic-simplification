#include "glue_face_sides.h"

void glue_face_sides(
    const Eigen::Vector2i & fs0,
    const Eigen::Vector2i & fs1,
    Eigen::MatrixXi & G)
    {
        assert(!is_same_face_side(fs0, fs1)); // cannot glue to itself
        using namespace global_variables;
        using namespace std;
        int f0 = fs0(0);
        int s0 = fs0(1);
        int f1 = fs1(0);
        int s1 = fs1(1);
        if (f0 == GHOST_INDEX && s0 == GHOST_INDEX && f1 != GHOST_INDEX && s1 != GHOST_INDEX)
        {
            // if fs1 is a boundary face side
            G.block(f1,s1*2,1,2) = fs0.transpose();
        }
        else if (f0 != GHOST_INDEX && s0 != GHOST_INDEX && f1 == GHOST_INDEX && s1 == GHOST_INDEX)
        {
            // if fs0 is a boundary face side
            G.block(f0,s0*2,1,2) = fs1.transpose();
        }
        else if (f0 != GHOST_INDEX && s0 != GHOST_INDEX && f1 != GHOST_INDEX && s1 != GHOST_INDEX)
        {
            // if interior edge
            G.block(f0,s0*2,1,2) = fs1.transpose();
            G.block(f1,s1*2,1,2) = fs0.transpose();
        }
        else
            assert(false && "gluing invalid face sides");
    }
#include "remove_ear_vertex.h"

int remove_ear_vertex(
    const int & v,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V,
    Eigen::MatrixXd & T)
{
    // TODO: this does not handle self-boundary-edge yet! MUST DO
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    assert (is_ear_vertex(G,v2fs,v));

    Vector2i fs = next(v2fs.row(v).transpose());
    Vector2i fs_twin = twin(G,fs);
    int f = fs(0);
    int f_twin = fs_twin(0);

    // check whether to-be-removed face still has barycentric points
    if (F2V[f].size() > 0)
    {
        cout << "[WARNING] ear face still has barycentric points. Mapping to its twin ...\n";

        // compute the 2D parameterization of the diamond
        MatrixXd U;
        MatrixXi UF;
		flatten_diamond_mesh(G,l,fs,U,UF);

        // get barycentric point indices
		vector<int> bIdx(F2V[f]); 

        // get barycentric points in 2D
		MatrixXd P(bIdx.size(), 2);
        {
			int pIdx = 0;
			RowVector2d u0, u1, u2;
			Vector3d bc;
			for (int ii=0; ii<F2V[f].size(); ii++)
			{
				int b = F2V[f][ii]; // bary index
				
				u0 << U.row(UF(0,0));
				u1 << U.row(UF(0,1));
				u2 << U.row(UF(0,2));
				bc << BC.row(b).transpose();
				
				P.row(pIdx) = bc(0)*u0 + bc(1)*u1 + bc(2)*u2;
				pIdx += 1;
			}
		}

        // compute barycentric coordiantes on f_twin
        for (int ii=0; ii<P.rows(); ii++)
		{
			// bary point
			Vector2d p;
			p << P.row(ii).transpose();

			// compute bary for on f_twin
			Vector2d u0, u1, u2;
			u0 << U.row(UF(1,0)).transpose();
			u1 << U.row(UF(1,1)).transpose();
			u2 << U.row(UF(1,2)).transpose();

			bool degen, in; // is degen & is in
			double d; // dist to valid
			Vector3d b = compute_barycentric_robust(p, u0, u1, u2, degen, in, d);

			// check whether is valud
            if (!(d < EPS)) // if outside 
				throw std::runtime_error("[Error] barycentric point is outside");

			// assign barycentric
			int bi = bIdx[ii];
            BC.row(bi) << b.transpose();
            F2V[f_twin].push_back(bi);
		}
        F2V[f].clear();
    }

    // store v as barycentric coordiantes in F2V[f_twin]
    MatrixXd U;
    MatrixXi UF;
    flatten_diamond_mesh(G,l,fs,U,UF);
    Vector2d u_v = U.row(2).transpose();

    bool degen, inside;
    double dist;
    Vector2d u0, u1, u2;
    u0 << U.row(UF(1,0)).transpose();
    u1 << U.row(UF(1,1)).transpose();
    u2 << U.row(UF(1,2)).transpose();
    Vector3d bary = compute_barycentric_robust(u_v, u0, u1, u2, degen, inside, dist);
    if (inside || (dist < EPS))
    {
        BC.row(v) << bary.transpose();
        F2V[f_twin].push_back(v);
    }
    else
        throw std::runtime_error("[Error] ear vertex has invalid barycentric coordinates, possibly due to having curvature");

    // update face sides to the valid one
    int v_right, v_left;
    get_face_side_vertices(F, fs, v_right, v_left);

    if (v_left != v_right) // if this is not a self boundary loop where v_left == v_right
    {
        // update v2fs for v_left
        v2fs.row(v_left) << fs_twin.transpose();

        // update v2fs for v_right by doing clockwise rotations until the end
        Vector2i fs_cw = fs;
        while (!is_same_face_side(twin(G,fs_cw), GHOST_FACE_SIDE)) 
            fs_cw = cw(G,fs_cw); // cw rotate to the right most face side
        v2fs.row(v_right) << fs_cw.transpose();
    }
    else
    {
        // this is a weird self boundary loop where the ear vertex look like (vi > vj > vi)
        // in this case, we assign v2fs to the self edge (whic is the fs_twin) after removal
        v2fs.row(v_left) << fs_twin.transpose();
    }

    // glue the twin of the fs opposite to v to boundary
    glue_face_sides(fs_twin, GHOST_FACE_SIDE, G);

    // update mesh info
    F.row(f) << GHOST_INDEX, GHOST_INDEX, GHOST_INDEX;
    l.row(f) << -1, -1, -1;
    G.row(f) << GHOST_INDEX, GHOST_INDEX, GHOST_INDEX, GHOST_INDEX, GHOST_INDEX, GHOST_INDEX;
    A.row(f) << DOUBLE_INF, DOUBLE_INF, DOUBLE_INF;

    return f_twin;
}

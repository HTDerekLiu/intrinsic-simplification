#include "flatten_ear_vertex_and_cost.h"

double pre_flatten_ear_vertex_and_cost(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    const Eigen::MatrixXd & T,
    const Eigen::VectorXd & K,
    const int & v_flatten)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    assert(is_ear_vertex(G,v2fs,v_flatten)); 

    vector<Vector2i> fs_list;
    vertex_one_ring_face_sides(G,v2fs,v_flatten,fs_list);
    assert(fs_list.size() == 1);

    // store histories in order to change them back
    Vector2i fs_flip = next(fs_list[0]);
    // cout << "fs_flip: " << fs_flip.transpose() << endl;

    // old edge lengths, each tuple stores (f, s, l_old)
    vector< tuple<int,int, double> > l_history; 
    // old face indices, each tuple stores (f, vi, vj, vk)
    vector< tuple<int, int,int,int> > F_history;
    // old glue map, each tuple stores (f,s,f_twin,s_twin)
    vector< tuple<int,int,int,int> > G_history;
    // old angular coordinates, each tuple stores (f, s, A_old)
    vector< tuple<int,int,double> > A_history;
    // old v2fs, each tuple stores (v, f, s)
    vector< tuple<int, int,int> > v2fs_history;
    {
        Vector2i fs_flip_twin = twin(G, fs_flip);
        if (is_same_face_side(fs_flip_twin, GHOST_FACE_SIDE)){
            // the opposite edge of the ear vertex is a boundary edge -> cannot be collapses
            return DOUBLE_INF;
        }

        int f = fs_flip(0);
        int s = fs_flip(1);
        int ft = fs_flip_twin(0);
        int st = fs_flip_twin(1);

        // add old edge lengths
        l_history.reserve(2 * 3);
        for (int ss=0; ss<3; ss++)
        {
            l_history.push_back(make_tuple(f, ss, l(f,ss)));
            l_history.push_back(make_tuple(ft, ss, l(ft,ss)));
        }
        
        // add old face indices
        F_history.reserve(2);
        F_history.push_back(make_tuple(f, F(f,0), F(f,1), F(f,2)));
        F_history.push_back(make_tuple(ft, F(ft,0), F(ft,1), F(ft,2)));

        // add glue map back
        G_history.reserve(3 * 2);
        for (int ss=0; ss<3; ss++)
        {
            Vector2i fss; fss << f, ss;
            Vector2i fss_twin; fss_twin << G(f,ss*2), G(f,ss*2+1);
            Vector2i ftss; ftss << ft, ss;
            Vector2i ftss_twin; ftss_twin << G(ft,ss*2), G(ft,ss*2+1);

            G_history.push_back(make_tuple(fss(0), fss(1), fss_twin(0), fss_twin(1)));
            G_history.push_back(make_tuple(ftss(0), ftss(1), ftss_twin(0), ftss_twin(1)));
        }

        // add angular coordinates
        A_history.reserve(2 * 3);
        for (int ss=0; ss<3; ss++)
        {
            A_history.push_back(make_tuple(f, ss, A(f,ss)));
            A_history.push_back(make_tuple(ft, ss, A(ft,ss)));
        }

        // add v2fs history
        int vi, vj, vk, vl, v_useless;
        get_face_side_vertices(F,fs_flip,vi,vj);
        get_face_side_vertices(F,next(fs_flip),v_useless,vk);
        get_face_side_vertices(F,next(twin(G,fs_flip)),v_useless,vl);
        v2fs_history.reserve(4);
        v2fs_history.push_back(make_tuple(vi, v2fs(vi,0), v2fs(vi,1)));
        v2fs_history.push_back(make_tuple(vj, v2fs(vj,0), v2fs(vj,1)));
        v2fs_history.push_back(make_tuple(vk, v2fs(vk,0), v2fs(vk,1)));
        v2fs_history.push_back(make_tuple(vl, v2fs(vl,0), v2fs(vl,1)));
    }

    // flip an edge (if it can be flipped to a doamond)
    if (is_diamond_convex(G,l,fs_flip))
        flip_edge(fs_flip, F,G,l,A,v2fs);
    else 
        return DOUBLE_INF;

    // flatten the boundary vertex
    double cost = pre_flatten_boundary_vertex_and_cost(F,G,l,A,v2fs,T,K,v_flatten);

    // put back the intrinsic info
    for (int ii=0; ii<l_history.size(); ii++)
    {
        auto fsl = l_history[ii];
        int f = get<0>(fsl);
        int s = get<1>(fsl);
        double l_fs = get<2>(fsl);
        l(f,s) = l_fs;
    }
    for (int ii=0; ii<F_history.size(); ii++)
    {
        auto fvivjvk = F_history[ii];
        int f = get<0>(fvivjvk);
        int vi = get<1>(fvivjvk);
        int vj = get<2>(fvivjvk);
        int vk = get<3>(fvivjvk);
        F.row(f) << vi, vj, vk;
    }
    for (int ii=0; ii<G_history.size(); ii++)
    {
        auto fsftst = G_history[ii];
        int f = get<0>(fsftst);
        int s = get<1>(fsftst);
        int ft = get<2>(fsftst);
        int st = get<3>(fsftst);
        Vector2i fs_tmp; fs_tmp << f, s;
        Vector2i fs_tmp_twin; fs_tmp_twin << ft, st;
        glue_face_sides(fs_tmp, fs_tmp_twin, G);
    }
    for (int ii=0; ii<A_history.size(); ii++)
    {
        auto fsA = A_history[ii];
        int f = get<0>(fsA);
        int s = get<1>(fsA);
        double angle = get<2>(fsA);
        A(f,s) = angle;
    }
    for (int ii=0; ii<v2fs_history.size(); ii++)
    {
        auto vfs = v2fs_history[ii];
        int v = get<0>(vfs);
        int f = get<1>(vfs);
        int s = get<2>(vfs);
        v2fs.row(v) << f, s;
    }
    return cost;
}

double flatten_ear_vertex_and_cost(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXi & v2fs,
    std::vector<std::vector<int>> & F2V,
    const int & v_flatten,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXd & BC,
    Eigen::MatrixXd & T,
    Eigen::VectorXd & K)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    assert(is_ear_vertex(G,v2fs,v_flatten)); 

    vector<Vector2i> fs_list;
    vertex_one_ring_face_sides(G,v2fs,v_flatten,fs_list);
    assert(fs_list.size() == 1);

    // store histories in order to change them back
    Vector2i fs_flip = next(fs_list[0]);

    // detect whether it can be flipped
    Vector2i fs_flip_twin = twin(G, fs_flip);
    if (is_same_face_side(fs_flip_twin, GHOST_FACE_SIDE)){
        // the opposite edge of the ear vertex is a boundary edge -> cannot be collapses
        return DOUBLE_INF;
    }

    // flip an edge (if it can be flipped to a doamond)
    if (is_diamond_convex(G,l,fs_flip))
        flip_edge(fs_flip, F,G,l,A,v2fs);
    else 
        return DOUBLE_INF;

    // flatten the boundary vertex
    double cost = flatten_boundary_vertex_and_cost(F,G,v2fs,F2V,v_flatten,l,A,BC,T,K);

    // flip edge back
    flip_edge_cw(fs_flip, F,G,l,A,v2fs,BC,F2V);

    return cost;
}

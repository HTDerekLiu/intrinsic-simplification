#include "flatten_interior_vertex_and_cost.h"

static void get_edge_lengths(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::MatrixXd & l_012,
	Eigen::MatrixXd & scale_012,
    std::vector<Eigen::Vector2i> & all_one_ring_fs)
{
	/*
	- each row of l_012 stores edge lengths from v0->v1, v1->v2, v2->v0 respectively
	- each row of scale_012 stores the multiple of scaling factor in l_012. In the usual case, each row of scale_012 will be [1, 0, 1] because only v0->v1 and v2->v0 touches v0. However, this won't be the case if this patch has self-edges
	*/
    using namespace std;
    using namespace Eigen;

    // vertex_one_ring_unique_face_sides(G, v2fs, v, one_ring_fs);
    vertex_one_ring_face_sides(G, v2fs, v, all_one_ring_fs);

    int nF_local = all_one_ring_fs.size();
    l_012.resize(nF_local, 3); 
	scale_012.resize(nF_local, 3);
	scale_012.setZero(); 
    for (int ii=0; ii<nF_local; ii++)
    {
		int vi, vj; 

        Vector2i fs = all_one_ring_fs[ii];
        l_012(ii, 0) = l(fs(0), fs(1));
		get_face_side_vertices(F,fs,vi,vj);
		if (vi == v)
			scale_012(ii, 0) += 1;
		if (vj == v)
			scale_012(ii, 0) += 1;
		
        Vector2i fs_next = next(fs);
        l_012(ii, 1) = l(fs_next(0), fs_next(1));
		get_face_side_vertices(F,fs_next,vi,vj);
		if (vi == v)
			scale_012(ii, 1) += 1;
		if (vj == v)
			scale_012(ii, 1) += 1;

        Vector2i fs_next_next = next(fs_next);
        l_012(ii, 2) = l(fs_next_next(0), fs_next_next(1));
		get_face_side_vertices(F,fs_next_next,vi,vj);
		if (vi == v)
			scale_012(ii, 2) += 1;
		if (vj == v)
			scale_012(ii, 2) += 1;
    }
}

double pre_flatten_interior_vertex_and_cost(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    const Eigen::MatrixXd & A,
    const Eigen::MatrixXi & v2fs,
    const Eigen::MatrixXd & T,
    const Eigen::VectorXd & K,
    const int & v_flatten)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    assert(is_interior_vertex(G,v2fs,v_flatten));
    assert(is_one_ring_fs_valid(F,G,v2fs,v_flatten));

    // gather edge lengths
    vector<Vector2i> one_ring_fs;
    MatrixXd l012, s012;
    get_edge_lengths(F, G,l,v2fs, v_flatten, l012, s012, one_ring_fs);

    // if interior vertices only have two one-ring fs 
    if (one_ring_fs.size() == 2)
        return DOUBLE_INF;

    // detect self faces
    for (int ii=0; ii<one_ring_fs.size(); ii++){
        int f = one_ring_fs[ii](0);
        if (F(f,0) == F(f,1) && F(f,0) == F(f,2)) // self face
            return DOUBLE_INF;
    }

    // CETM flattening 
    double u; // conformal scaling factor
    cetm_flatten_interior_vertex(l012, s012, u);
    if (isnan(u))
        return DOUBLE_INF;
    
    // once you go pass this point, the coarsening MUST not fail in the later procedure theoritically. For efficient implementation, let's change the mesh and then change it back

    // get one-ring vertices
    VectorXi one_ring_v;
    vector<Vector2i> fs_of_one_ring_v;
    bool _; 
    vertex_one_ring_vertices(F,G,v2fs, v_flatten, one_ring_v, _, fs_of_one_ring_v);

    // compute the Gaussian curvature before scaling
    int num_one_ring_v = one_ring_v.size();
    VectorXd Ks_pre(num_one_ring_v);
    {
        for (int ii=0; ii<num_one_ring_v; ii++)
        {
            int vj = one_ring_v(ii);
            Ks_pre(ii) = K(vj);
        }
    }

    // get one ring faces
    int num_one_ring_fs = one_ring_fs.size();
    unordered_set<int> one_ring_f_set;
    for (int ii=0; ii<num_one_ring_fs; ii++)
        one_ring_f_set.insert(one_ring_fs[ii](0));
    
    // store the unscaled edge lengths in order to change it back later
    int num_one_ring_f = one_ring_f_set.size();
    vector< pair<int, RowVector3d> > l_history(num_one_ring_f); 
    {
        int count = 0;
        for (const int& f : one_ring_f_set)
        {
            l_history[count] = make_pair(f,l.row(f));
            count ++;
        }
    }

    // scale the edge lengths
    for (const int& f : one_ring_f_set)
    {
		int vi, vj; 
        for (int s=0; s<3; s++){
            Vector2i fs; 
            fs << f, s;
            get_face_side_vertices(F,fs,vi,vj);
            
            double multiple = 0.0;
            if (vi == v_flatten)
                multiple += 1.0;
            if (vj == v_flatten)
                multiple += 1.0;
            l(f, s) *= exp(multiple * u / 2.0);
        }
    }
    assert(abs(gaussian_curvature_at_vertex(G,l,v2fs,v_flatten)) < 1e-7);

    // compute the Guassian curvature after scaling
    VectorXd Ks_post(num_one_ring_v);
    {
        for (int ii=0; ii<num_one_ring_v; ii++)
        {
            int vj = one_ring_v(ii);
            Ks_post(ii) = gaussian_curvature_at_vertex(G, l, v2fs, vj);
        }
    }

    // compute the transport ratio
    VectorXd dK = (Ks_post - Ks_pre).cwiseAbs();
    if (dK.sum() == 0.0) // avoid all zeros
        dK = dK.array() + 1e-6;
    
    // compute transport cost and update transport info T
    int num_channels = T.cols() / 3; 
    double cost = 0.0;
    {
        for (int ii=0; ii<fs_of_one_ring_v.size(); ii++)
        {
            Vector2i fs;
            fs << fs_of_one_ring_v[ii];
            double alpha_ij = dK(ii) / dK.sum();
            int vj = one_ring_v[ii];

            for (int ch=0; ch<num_channels; ch++)
            {
                // compute curvature transport at the tangent plane of vi
                // warning: T is using the global index
                double alpha_mi = T(v_flatten, ch*3) * alpha_ij;
                Vector2d ti; ti << T(v_flatten, ch*3+1), T(v_flatten, ch*3+2);

                // if nothing to transfer, then skip this for loop
                if (alpha_mi == 0) 
                    continue;

                // parallel transport to the tangent plane of vj (only rotation)
                Vector2i fs_twin = twin(G, fs);
                double theta = A(fs_twin(0), fs_twin(1)) - A(fs(0), fs(1)) + M_PI;
                double c = cos(theta);
                double s = sin(theta);
                MatrixXd Rij(2, 2);
                Rij << c, -s,
                       s,  c;

                // compute the edge vector eji = -eij at the tangent plane of vj
                Vector2d logj_eji;
                double l_fs_twin = l(fs_twin(0), fs_twin(1));
                double A_fs_twin = A(fs_twin(0), fs_twin(1));
                polar_to_cartesian(l_fs_twin, A_fs_twin, logj_eji);

                // combine transported curvature with existing curvature at vj
                // Note: since this is preflattening, we don't change T
                double mj = T(vj, ch*3);
                Vector2d tj; tj << T(vj, ch*3+1), T(vj, ch*3+2);
                
                // compute cost
                double mj_new = alpha_mi + mj;
                Vector2d tj_new; 
                tj_new = (alpha_mi*(Rij*ti + logj_eji) + mj*tj) / (alpha_mi + mj);
                cost += mj_new * tj_new.norm();
            }
        }
    }

    // scale the edge lengths back to the original
    for (int ii=0; ii<l_history.size(); ii++)
    {
        auto fl = l_history[ii];
        int f = fl.first;
        RowVector3d l_row_f = fl.second;
        l.row(f) << l_row_f;
    }

    return cost;
}


double flatten_interior_vertex_and_cost(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const std::vector<std::vector<int>> & F2V,
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

    assert (is_interior_vertex(G,v2fs,v_flatten));
    assert(is_one_ring_fs_valid(F,G,v2fs,v_flatten));

    // gather edge lengths
    vector<Vector2i> one_ring_fs;
    MatrixXd l012, s012;
    get_edge_lengths(F, G,l,v2fs, v_flatten, l012, s012, one_ring_fs);

    // if interior vertices only have two one-ring fs 
    if (one_ring_fs.size() == 2)
        return DOUBLE_INF;

    // detect self faces
    for (int ii=0; ii<one_ring_fs.size(); ii++){
        int f = one_ring_fs[ii](0);
        if (F(f,0) == F(f,1) && F(f,0) == F(f,2)) // self face
            return DOUBLE_INF;
    }

    // CETM flattening 
    double u; // conformal scaling factor
    cetm_flatten_interior_vertex(l012, s012, u);
    if (isnan(u))
        return DOUBLE_INF;
    
    // once you go pass this point, the coarsening MUST not fail in the later procedure theoriticall, so let's commit and change the output mesh

    // get one-ring vertices
    VectorXi one_ring_v;
    vector<Vector2i> fs_of_one_ring_v;
    bool _; 
    vertex_one_ring_vertices(F,G,v2fs, v_flatten, one_ring_v, _, fs_of_one_ring_v);

    // compute the Gaussian curvature before scaling
    int num_one_ring_v = one_ring_v.size();
    VectorXd Ks_pre(num_one_ring_v);
    {
        for (int ii=0; ii<num_one_ring_v; ii++)
        {
            int vj = one_ring_v(ii);
            Ks_pre(ii) = K(vj);
        }
    }

    // get one ring faces
    int num_one_ring_fs = one_ring_fs.size();
    unordered_set<int> one_ring_f_set;
    for (int ii=0; ii<num_one_ring_fs; ii++)
        one_ring_f_set.insert(one_ring_fs[ii](0));

    // scale the edge lengths
    for (const int& f : one_ring_f_set)
    {
		int vi, vj; 
        for (int s=0; s<3; s++){
            Vector2i fs; 
            fs << f, s;
            get_face_side_vertices(F,fs,vi,vj);
            
            double multiple = 0.0;
            if (vi == v_flatten)
                multiple += 1.0;
            if (vj == v_flatten)
                multiple += 1.0;
            l(f, s) *= exp(multiple * u / 2.0);
        }
    }
    assert(abs(gaussian_curvature_at_vertex(G,l,v2fs,v_flatten)) < 1e-7);

    // update angular coordinates (v_flatten's one ring)
    for (int jj=0; jj<num_one_ring_v; jj++)
    {
        int vj = one_ring_v(jj);
        update_angular_coordinate(G,l,v2fs,vj,A);
    }
     
    // compute the Guassian curvature after scaling
    VectorXd Ks_post(num_one_ring_v);
    {
        for (int ii=0; ii<num_one_ring_v; ii++)
        {
            int vj = one_ring_v(ii);
            double Kvj = gaussian_curvature_at_vertex(G, l, v2fs, vj);
            Ks_post(ii) = Kvj;
            K(vj) = Kvj;
        }
    }

    // compute the transport ratio
    VectorXd dK = (Ks_post - Ks_pre).cwiseAbs();
    if (dK.sum() == 0.0) // avoid all zeros
        dK = dK.array() + 1e-6;

    // update the barycnetric via projective interpolation
    for (const int& ff : one_ring_f_set) // loop one-ring faces
    {
        for (int jj=0; jj<F2V[ff].size(); jj++) // get index to BC
        {
            int bIdx = F2V[ff][jj];
            for (int s=0; s<3; s++)
            {
                if (F(ff,s) == v_flatten)
                {
                    BC(bIdx,s) *= exp(u); // projectively scale the barycentric
                }
            }
            BC.row(bIdx) /= BC.row(bIdx).sum(); 
        }
    }

    // compute transport cost and update transport info T
    int num_channels = T.cols() / 3; 
    double cost = 0.0;
    {
        for (int ii=0; ii<fs_of_one_ring_v.size(); ii++)
        {
            Vector2i fs;
            fs << fs_of_one_ring_v[ii];
            double alpha_ij = dK(ii) / dK.sum();
            int vj = one_ring_v[ii];

            // get the opposite vertex vj
            int useless, vj_test;
            get_face_side_vertices(F, fs, useless, vj_test);

            for (int ch=0; ch<num_channels; ch++)
            {
                // compute curvature transport at the tangent plane of vi
                // warning: T is using the global index
                double alpha_mi = T(v_flatten, ch*3) * alpha_ij;
                Vector2d ti; ti << T(v_flatten, ch*3+1), T(v_flatten, ch*3+2);

                // if nothing to transfer, then skip this for loop
                if (alpha_mi == 0) 
                    continue;

                // parallel transport to the tangent plane of vj (only rotation)
                Vector2i fs_twin = twin(G, fs);
                double theta = A(fs_twin(0), fs_twin(1)) - A(fs(0), fs(1)) + M_PI;
                double c = cos(theta);
                double s = sin(theta);
                MatrixXd Rij(2, 2);
                Rij << c, -s,
                       s,  c;

                // compute the edge vector eji = -eij at the tangent plane of vj
                Vector2d logj_eji;
                double l_fs_twin = l(fs_twin(0), fs_twin(1));
                double A_fs_twin = A(fs_twin(0), fs_twin(1));
                polar_to_cartesian(l_fs_twin, A_fs_twin, logj_eji);

                // combine transported curvature with existing curvature at vj
                // Note: since this is preflattening, we don't change T
                double mj = T(vj, ch*3);
                Vector2d tj; tj << T(vj, ch*3+1), T(vj, ch*3+2);
                
                // compute cost
                double mj_new = alpha_mi + mj;
                Vector2d tj_new; 
                tj_new = (alpha_mi*(Rij*ti + logj_eji) + mj*tj) / (alpha_mi + mj);
                cost += mj_new * tj_new.norm();

                // update T
                T(vj, ch*3) = mj_new;
                T(vj, ch*3+1) = tj_new(0);
                T(vj, ch*3+2) = tj_new(1);
            }
        }
        T.row(v_flatten).setZero();
    }
    return cost;
}
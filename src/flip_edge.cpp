#include "flip_edge.h"

static void relable(
	const Eigen::Vector2i & s0,
	const Eigen::Vector2i & s1,
	const Eigen::Vector2i & s2,
	const Eigen::Vector2i & t0,
	const Eigen::Vector2i & t1,
	const Eigen::Vector2i & t2,
	Eigen::Vector2i & s)
{
	// relabel the face side to handle self edges. This relabeling is required because the twin of a boundary face side will be inside the diamond. And the face side inside the diamond may get changed. So we have to make sure it is correct via relabeling!
	using namespace std;
	if (is_same_face_side(s, s1))
	{
		s = next(next(t0));
		return;
	}
	if (is_same_face_side(s, s2))
	{
		s = next(s0);
		return;
	}
	if (is_same_face_side(s, t1))
	{
		s = next(next(s0));
		return;
	}
	if (is_same_face_side(s, t2))
	{
		s = next(t0);
		return;
	}
}

void flip_edge(
    const Eigen::Vector2i & s0,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V,
    std::map<std::pair<int, int>, std::pair<int, int>> & fs_dict)
{
    using namespace Eigen;
	using namespace std;
	using namespace global_variables;

	assert(not is_boundary_face_side(G,s0));
	assert(is_diamond_convex(G,l,s0));

	Vector2i t0 = twin(G, s0);
	Vector2i t1 = next(t0);
	Vector2i t2 = next(t1);

	Vector2i s1 = next(s0);
	Vector2i s2 = next(s1);
	
	Vector2i s1t = twin(G, s1);
	Vector2i s2t = twin(G, s2);
	Vector2i t1t = twin(G, t1);
	Vector2i t2t = twin(G, t2);

	double angle_s1 = get_face_side_angular_coordinate(A, s1);
	double angle_s2 = get_face_side_angular_coordinate(A, s2);
	double angle_t1 = get_face_side_angular_coordinate(A, t1);
	double angle_t2 = get_face_side_angular_coordinate(A, t2);

    // create a list of face side keys (used to create a map later)
	Vector2i s0_key, s1_key, s2_key, t0_key, t1_key, t2_key, s1t_key, s2t_key, t1t_key, t2t_key;
	{
		s0_key << s0;
		s1_key << s1;
		s2_key << s2;
		t0_key << t0;
		t1_key << t1;
		t2_key << t2;
		s1t_key << s1t;
		s2t_key << s2t;
		t1t_key << t1t;
		t2t_key << t2t;
	}

    // Get vertex indices for the vertices of the diamond
	int v0, v1, v2, v3;
	{
		v0 = F(s0(0), s0(1));
		v1 = F(s1(0), s1(1));
		v2 = F(s2(0), s2(1));
		v3 = F(t2(0), t2(1));
	}

    // Get the two faces from our face-sides
    int f0, f1;
	{
		f0 = s0(0);
		f1 = t0(0);
	}

    // Get the original lengths of the outside edges of the diamond
	double ls1, ls2, lt1, lt2;
	{
		ls1 = l(s1(0), s1(1));
		ls2 = l(s2(0), s2(1));
		lt1 = l(t1(0), t1(1));
		lt2 = l(t2(0), t2(1));
	}

    // Compute the length of the new edge
    double new_length = diagonal_length(G, l, s0);

    // assign the new vertex indices and edge lengths
	Vector3i f0_vIdx, f1_vIdx;
	{
		Vector3i tmp; 
		tmp << v3, v2, v0;
		roll1d(tmp, s0(1), f0_vIdx);
		tmp << v2, v3, v1;
		roll1d(tmp, t0(1), f1_vIdx);
	}
	Vector3d f0_l, f1_l;
	{
		Vector3d tmp; 	
		tmp << new_length, ls2, lt1;
		roll1d(tmp, s0(1), f0_l);
		tmp << new_length, lt2, ls1;
		roll1d(tmp, t0(1), f1_l);
	}

    // update barycentric coordinates   
	if (F2V[f0].size()>0 || F2V[f1].size()>0) // if there are barycentric points in this diamond
	{ 
		// compute the 2D parameterization of the diamond
        MatrixXd U;
        MatrixXi UF,UFF;
		flatten_diamond_mesh(G,l,s0,true,U,UF,UFF);
		// if (f0 == 11 && f1 == 53){
		// 	cout << "U: \n" << U << endl;
		// 	cout << "UF: \n" << UF << endl;
		// 	cout << "UFF: \n" << UFF << endl;
		// }
	
		// get barycentric point indices
		vector<int> bIdx(F2V[f0]); // combined bary indices
		{
    		bIdx.insert(bIdx.end(), F2V[f1].begin(), F2V[f1].end());
		}
		
		// get barycentric points in 2D
		MatrixXd P(bIdx.size(), 2);
		{
			int pIdx = 0;
			RowVector2d u0, u1, u2;
			Vector3d bc;
			for (int ii=0; ii<F2V[f0].size(); ii++)
			{
				int b = F2V[f0][ii]; // bary index
				
				u0 << U.row(UF(0,0));
				u1 << U.row(UF(0,1));
				u2 << U.row(UF(0,2));
				bc << BC.row(b).transpose();
				
				P.row(pIdx) = bc(0)*u0 + bc(1)*u1 + bc(2)*u2;
				pIdx += 1;
			}
			for (int ii=0; ii<F2V[f1].size(); ii++)
			{
				int b = F2V[f1][ii]; // bary index

				u0 << U.row(UF(1,0));
				u1 << U.row(UF(1,1));
				u2 << U.row(UF(1,2));
				bc << BC.row(b).transpose();
				
				P.row(pIdx) = bc(0)*u0 + bc(1)*u1 + bc(2)*u2;
				pIdx += 1;
			}
			assert(!isnan(bc(0)));
			assert(!isnan(bc(1)));
			assert(!isnan(bc(2)));
		}

		// compute barycentric coordinate of the new mesh
		vector<int> bIdx_f0, bIdx_f1;

		// reserve bIdx.size, that is the max possible size
		bIdx_f0.reserve(bIdx.size()); 
		bIdx_f1.reserve(bIdx.size());

		for (int ii=0; ii<P.rows(); ii++)
		{
			// bary point
			Vector2d p;
			p << P.row(ii).transpose();

			// compute bary for UFF0
			Vector2d u0_0, u0_1, u0_2;
			u0_0 << U.row(UFF(0,0)).transpose();
			u0_1 << U.row(UFF(0,1)).transpose();
			u0_2 << U.row(UFF(0,2)).transpose();

			bool degen0, in0; // is degen & is in
			double d0; // dist to valid
			Vector3d b0 = compute_barycentric_robust(p, u0_0, u0_1, u0_2, degen0, in0, d0);

			// compute bary for UFF1
			Vector2d u1_0, u1_1, u1_2;
			u1_0 << U.row(UFF(1,0)).transpose();
			u1_1 << U.row(UFF(1,1)).transpose();
			u1_2 << U.row(UFF(1,2)).transpose();

			bool degen1, in1;
			double d1;
			Vector3d b1 = compute_barycentric_robust(p, u1_0, u1_1, u1_2, degen1, in1, d1);

			// determin f0 or f1
			bool to_f0;
			if (!degen0 && !degen1) // f0,f1 not degenerated
			{
				if (in0 && !in1) // inside f0
					to_f0 = true;
				else if (!in0 && in1) // inside f1
					to_f0 = false;
				else if (in0 && in1) // on the edge
					if (d0 <= d1) // closer to f0
						to_f0 = true;
					else // closer to f1
						to_f0 = false;
				else // outside the mesh
				{
					cout << "barycentric points on face " << f0 << " or " << f1 << " is outside\n";
					cout << "b0: " << b0.transpose() << endl;
					cout << "b1: " << b1.transpose() << endl;
					throw std::runtime_error("[Error] barycentric point is outside (no degenerated faces)");
				}
			}
			else if (!degen0 && degen1) // only f1 degenerated
			{
				if (in0 || (d1<EPS)) // inside f0
					to_f0 = true;
				else // not inside f0
					throw std::runtime_error("[Error] barycentric point is outside (f1 degenerated)");
			}
			else if (degen0 && !degen1) // f0 is degenerated
			{
				if (in1 || (d1<EPS)) // inside f1
					to_f0 = false;
				else // not in f1
					throw std::runtime_error("[Error] barycentric point is outside (f0 degenerated)");
			}
			else
			{
				cout << "U: \n" << U << endl;
				cout << "UFF: \n" << UFF << endl;
				throw std::runtime_error("[Error] entire diamond is degenerated");
			}

			// assign barycentric
			int bi = bIdx[ii];
			if (to_f0)
			{
				bIdx_f0.push_back(bi);
				BC.row(bi) << b0.transpose();
			}
			else
			{
				bIdx_f1.push_back(bi);
				BC.row(bi) << b1.transpose();
			}
		}

		// modify F2V
		F2V[f0] = bIdx_f0;
		F2V[f1] = bIdx_f1;
	} 

    // Update the face list
    F.row(f0) = f0_vIdx;
    F.row(f1) = f1_vIdx;

	// Update the gluing map G
	relable(s0, s1, s2, t0, t1, t2, s1t);
	relable(s0, s1, s2, t0, t1, t2, s2t);
	relable(s0, s1, s2, t0, t1, t2, t1t);
	relable(s0, s1, s2, t0, t1, t2, t2t);

    // create a map from old face side to new face side
	Vector2i s1_new, s2_new, t1_new, t2_new;
	{
		s1_new = next(next(t0));
		s2_new = next(s0);
		t1_new = next(next(s0));
		t2_new = next(t0);
	}
    fs_dict.clear();
	fs_dict.insert({make_pair(s0_key(0), s0_key(1)), make_pair(s0(0), s0(1))});
	fs_dict.insert({make_pair(t0_key(0), t0_key(1)), make_pair(t0(0), t0(1))});
	fs_dict.insert({make_pair(s1_key(0), s1_key(1)), make_pair(s1_new(0), s1_new(1))});
	fs_dict.insert({make_pair(s2_key(0), s2_key(1)), make_pair(s2_new(0), s2_new(1))});
	fs_dict.insert({make_pair(t1_key(0), t1_key(1)), make_pair(t1_new(0), t1_new(1))});
	fs_dict.insert({make_pair(t2_key(0), t2_key(1)), make_pair(t2_new(0), t2_new(1))});
	fs_dict.insert({make_pair(s1t_key(0), s1t_key(1)), make_pair(s1t(0), s1t(1))});
	fs_dict.insert({make_pair(s2t_key(0), s2t_key(1)), make_pair(s2t(0), s2t(1))});
	fs_dict.insert({make_pair(t1t_key(0), t1t_key(1)), make_pair(t1t(0), t1t(1))});
	fs_dict.insert({make_pair(t2t_key(0), t2t_key(1)), make_pair(t2t(0), t2t(1))});

    glue_face_sides(s1_new, s1t, G);
    glue_face_sides(s2_new, s2t, G);
    glue_face_sides(t1_new, t1t, G);
    glue_face_sides(t2_new, t2t, G);

    // assign the face side angular coordinates (except (f0,0) and (f1,0))
	{
		if (!is_same_face_side(s2_new, GHOST_FACE_SIDE))
			A(s2_new(0), s2_new(1)) = angle_s2;
		if (!is_same_face_side(t1_new, GHOST_FACE_SIDE))
			A(t1_new(0), t1_new(1)) = angle_t1;
		if (!is_same_face_side(t2_new, GHOST_FACE_SIDE))
			A(t2_new(0), t2_new(1)) = angle_t2;
		if (!is_same_face_side(s1_new, GHOST_FACE_SIDE))
			A(s1_new(0), s1_new(1)) = angle_s1;
	}


    // Update the edge length because even the edges we didn't flip, they have been re-labeled, so we need to update those too.
    l.row(f0) = f0_l;
    l.row(f1) = f1_l;

	// update angular coordinate
	{
		Vector2i s0_twin = twin(G, s0);
		insert_angular_coordinate(G,l,s0,A);
		insert_angular_coordinate(G,l,s0_twin,A);
	}

    // reassign v2fs, in case the reference face-side of a vertex gets changed 
	MatrixXi fs_to_update(4,2); // each row is a fs to be checked
	{
		Vector2i fs = s0;
		Vector2i fs_next = next(fs);
		Vector2i fs_next_next = next(fs_next);
		Vector2i fs_twin_next_next = next(next(twin(G, fs)));

		fs_to_update.row(0) << fs(0), fs(1);
		fs_to_update.row(1) << fs_next(0), fs_next(1);
		fs_to_update.row(2) << fs_next_next(0), fs_next_next(1);
		fs_to_update.row(3) << fs_twin_next_next(0), fs_twin_next_next(1);
	}
	std::vector<int> vIdx;
	for (int ii=0; ii<fs_to_update.rows(); ii++)
	{
		Vector2i fs_cur = fs_to_update.row(ii).transpose();
		if (!is_same_face_side(fs_cur, GHOST_FACE_SIDE)) // if not outside
		{
			int vi = F(fs_cur(0), fs_cur(1));
			Vector2i fs_min = get_smallest_angular_coordinate(G, A, fs_cur);
			v2fs.row(vi) = fs_min.transpose();
			vIdx.push_back(vi);
		}
	}
}

void flip_edge(
    const Eigen::Vector2i & s0,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V)
{
	std::map<std::pair<int, int>, std::pair<int, int>> fs_dict;
	flip_edge(s0,F,G,l,A,v2fs,BC,F2V,fs_dict);
}

void flip_edge(
    const Eigen::Vector2i & s0,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs)
{
	std::map<std::pair<int, int>, std::pair<int, int>> fs_dict;
	Eigen::MatrixXd BC;
	std::vector<std::vector<int>> F2V(F.rows());
	flip_edge(s0,F,G,l,A,v2fs,BC,F2V,fs_dict);
}
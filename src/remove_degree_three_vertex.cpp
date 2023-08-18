#include "remove_degree_three_vertex.h"

static void relable(
	const std::vector<Eigen::Vector2i> & boundary_fs,
    const int & f_new,
	std::vector<Eigen::Vector2i> & boundary_fs_twin)
{
	using namespace std;
    using namespace Eigen;
    assert(boundary_fs.size() == 3);
    assert(boundary_fs_twin.size() == 3);

    Vector2i bfs0 = boundary_fs[0];
    Vector2i bfs1 = boundary_fs[1];
    Vector2i bfs2 = boundary_fs[2];
    for (int ii=0; ii<3; ii++)
    {
        if (is_same_face_side(bfs0, boundary_fs_twin[ii])){
            // cout << "relabel: "  << boundary_fs_twin[ii].transpose() << f_new << "  0" << endl;
            boundary_fs_twin[ii] << f_new, 0;
            continue;
        }
        else if (is_same_face_side(bfs1, boundary_fs_twin[ii])){
            // cout << "relabel: "  << boundary_fs_twin[ii].transpose() << f_new << "  1" << endl;
            boundary_fs_twin[ii] << f_new, 1;
            continue;
        }
        else if (is_same_face_side(bfs2, boundary_fs_twin[ii])){
            // cout << "relabel: "  << boundary_fs_twin[ii].transpose() << f_new << "  2" << endl;
            boundary_fs_twin[ii] << f_new, 2;
            continue;
        }
    }
}

int remove_degree_three_vertex(
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
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    assert (is_interior_vertex(G,v2fs,v));
    assert(is_one_ring_fs_valid(F,G,v2fs,v));

    vector<Vector2i> one_ring_fs;
    vertex_one_ring_face_sides(G, v2fs, v, one_ring_fs);
    assert (one_ring_fs.size() == 3);

    // gather local information all at once
    vector<int> faces(3); // 3 neighboring faces
    vector<int> vertices(3); // 3 one-ring vertices
    vector<double> bd_edge_lengths(3); // three opposite edge lengths
    vector<Vector2i> boundary_fs(3); // opposite face sides (v1,v2), (v2,v3), (v3,v1)
    vector<Vector2i> boundary_fs_twin(3); // twin of boundary_fs
    vector<double> angles_boundary_fs(3); // corner angles of boundary_fs
    for (int ii=0; ii<3; ii++)
    {
        Vector2i fs = one_ring_fs[ii];
        Vector2i fs_twin = twin(G,fs);
        Vector2i fs_next = next(fs);

        faces[ii] = fs(0);
        vertices[ii] = F(fs_twin(0), fs_twin(1));
        // vertices[ii] = F(fs_next(0), fs_next(1));
        bd_edge_lengths[ii] = l(fs_next(0), fs_next(1));
        boundary_fs[ii] = fs_next;
        angles_boundary_fs[ii] = A(fs_next(0), fs_next(1));
        boundary_fs_twin[ii] = twin(G, fs_next);
    }

    // compute the barycentric coodinate of the degree three vertex
    Vector3i bv;
    Vector3d bc; 
    {
        double area_sum = 0.0;
        double area;
        for (int ii=0; ii<3; ii++)
        {
            Vector2i fs = one_ring_fs[ii];
            Vector2i fs_next = next(fs);
            Vector2i fs_next_next = next(fs_next);

            // extract edge lengths
            double e0 = l(fs(0), fs(1));
            double e1 = l(fs_next(0), fs_next(1));
            double e2 = l(fs_next_next(0), fs_next_next(1));

            // compute face areas
            double s = (e0+e1+e2) / 2.0;
            if ((s*(s-e0)*(s-e1)*(s-e2)) < 0)
                area = 0.0;
            else
                area = sqrt(s * (s - e0) * (s - e1) * (s - e2));
            
            // keep track area sum 
            area_sum += area;

            // save barycentric info
            bc(ii) = area;
            Vector2i fs_opp = twin(G,next(next(twin(G,next(next(fs)))))); // the opposite vertex of this face 
            int _, vj;
            get_face_side_vertices(F, fs_opp, _, vj);
            bv(ii) = vj;
        }
        
        if (area_sum > 0)
            bc = bc / area_sum; // if triangle not degenerated, normalize bc
        else{
            // if triangle is degenerated, assign the barycentric to one of the edges (default edge 01)
            area_sum = 0; 
            Vector2i fs0 = one_ring_fs[0];
            Vector2i fs1 = one_ring_fs[1];
            Vector2i fs2 = one_ring_fs[2];

            double e0 = l(fs0(0), fs0(1));
            double e1 = l(fs1(0), fs1(1));

            bc(0) = e1 / (e0 + e1);
            bc(1) = e0 / (e0 + e1);
            bc(2) = 0.0;

            int v0, v1, v2, _;
            get_face_side_vertices(F, fs0, _, v0);
            get_face_side_vertices(F, fs1, _, v1);
            get_face_side_vertices(F, fs2, _, v2);
            bv << v0, v1, v2;
        }

        // shift bc bv
        int shift;
        for (int ii=0; ii<3; ii++)
            if ( vertices[ii] == bv(0) )
                shift = ii;
        Vector3d bc_unshift = bc;
        roll1d(bc_unshift, shift, bc);
        Vector3i bv_unshift = bv;
        roll1d(bv_unshift, shift, bv);

        // put back to BC
        BC.row(v) << bc.transpose();        
    }

    // update the barycentric coordinates of all existing barycentric coordinates on the removed faces
    vector<int> to_append_v_to_f;
    to_append_v_to_f.emplace_back(v);
    for (int i=0; i<faces.size(); i++)
    {
        int f = faces[i];
        for (int j=0; j<F2V[f].size(); j++)
        {
            int b = F2V[f][j];
            RowVector3d bc_b; bc_b.setZero();
            for (int ii=0; ii<3; ii++)
            {
                double bary_coord = BC(b,ii);
                int bary_vert = F(f,ii);
                if (bary_vert == v)
                    bc_b += bary_coord * bc;
                else
                {
                    int idx_in_bv = -1;
                    for (int k=0; k<3; k++)
                        if ( bv(k) == bary_vert )
                            idx_in_bv = k;
                    bc_b(idx_in_bv) += bary_coord;
                }
            }
            BC.row(b) = bc_b / bc_b.sum();
            to_append_v_to_f.emplace_back(b);
        } 
        F2V[f].clear();
    }
    // cout << "done bary update\n";

    // the face index of the new face
    int f_new = *std::min_element(faces.begin(), faces.end());
    // cout << "done finding new fIdx: " << f_new << "\n";

    // remove deleted face info
    for (int ii=0; ii<faces.size(); ii++)
    {
        int f = faces[ii];
        F.row(f) << GHOST_INDEX, GHOST_INDEX, GHOST_INDEX;
        l.row(f) << -1, -1, -1;
        G.row(f) << GHOST_INDEX, GHOST_INDEX, GHOST_INDEX, GHOST_INDEX, GHOST_INDEX, GHOST_INDEX;
        A.row(f) << DOUBLE_INF, DOUBLE_INF, DOUBLE_INF;
    }
    v2fs.row(v) << GHOST_INDEX, GHOST_INDEX;
    // cout << "done remove info\n";
    
    // add new face
    F.row(f_new) << vertices[0], vertices[1], vertices[2];

    // relabel boundary twin (so that self edges can be handled properly)
    relable(boundary_fs, f_new, boundary_fs_twin);

    // assign a barycentric coordinate the new face
    F2V[f_new] = to_append_v_to_f;

    // update angylar coordinates
    for (int ii=0; ii<3; ii++)
        A(f_new, ii) = angles_boundary_fs[ii];

    // glue the boundary_fs_twins with new face sides
    for (int ii=0; ii<3; ii++)
    {
        Vector2i fs; 
        fs << f_new, ii;
        glue_face_sides(boundary_fs_twin[ii], fs, G);
    }
    // cout << "done glue\n";

    // update edge length lists
    for (int ii=0; ii<3; ii++)
        l(f_new, ii) = bd_edge_lengths[ii];

    // update v2fs so that v2fs[v] returns the fs with smallest angular coordinate
    for (int ii=0; ii<3; ii++)
    {
        Vector2i fs; 
        fs << f_new, ii;
        Vector2i fs_min = get_smallest_angular_coordinate(G, A, fs);

        int vv = F(f_new, ii);
        v2fs.row(vv) << fs_min(0), fs_min(1);
    }

    return f_new;
}

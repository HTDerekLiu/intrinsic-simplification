#include "vertex_one_ring_vertices.h"

void vertex_one_ring_vertices(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::VectorXi & one_ring_vertices,
    bool & is_boundary_vertex,
    std::vector<Eigen::Vector2i> & face_sides_of_one_ring_vertices)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    is_boundary_vertex = false;
    vector<Vector2i> fs_list;
    vertex_one_ring_face_sides(G, v2fs, v, fs_list);

    if (fs_list.size() == 0) // if vertex v does not exist
    {
        one_ring_vertices = VectorXi();
        return;
    }

    // initialize a unordered set to keep track of unique vertices
	unordered_set<int> unique_vertex_set;
	unique_vertex_set.insert(v); // put the center vertex in, so that it won't get inserted

    // reserve a vertex list
    int nF = fs_list.size();
    vector<int> v_list;
    v_list.reserve(nF + 1); // +1 is for boundary vertices
    face_sides_of_one_ring_vertices.clear();
    face_sides_of_one_ring_vertices.reserve(nF + 1);
    
    int vi, vj;
    for (int ii=0; ii<nF; ii++)
    {
        Vector2i fs = fs_list[ii];
        get_face_side_vertices(F, fs, vi, vj);

        if (unique_vertex_set.find(vj) == unique_vertex_set.end())
        {
            // if vj does not exist
            v_list.push_back(vj);
            face_sides_of_one_ring_vertices.push_back(fs);
            unique_vertex_set.insert(vj);
        }

        // if encounter boundary face side
        if (is_same_face_side(ccw(G,fs), GHOST_FACE_SIDE))
        {
            Vector2i bd_fs_twin = next(next(fs)); // we can only use its twin because bd_fs is outside the model
            get_face_side_vertices(F, bd_fs_twin, vi, vj);

            // Note:
            // it is tempting to check whether vi has been added here like
            // >> if (unique_vertex_set.find(vi) == unique_vertex_set.end()) 
            // However, if we add this condition, the face_sides_of_one_ring_vertices may miss a face side if this is a "two-vertex boundary loop" which are two half-edges from vj -> vi and from vj -> vi. They only have two vertices, but they form a boundary loop. Thus we choose to allow this function to report a redundant vi for two-vertex boundary loop. This will make our transport cost computation a bit in accurate (because we count vi twice), but this can make sure our connectivity maintains valid.
            v_list.push_back(vi);
            face_sides_of_one_ring_vertices.push_back(bd_fs_twin); 
            unique_vertex_set.insert(vi);
            
            is_boundary_vertex = true;
        }
    }

    // put vertices from set into a vector
    one_ring_vertices.resize(v_list.size());
    one_ring_vertices = VectorXi::Map(v_list.data(), v_list.size());
}

void vertex_one_ring_vertices(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::VectorXi & one_ring_vertices,
    bool & is_boundary_vertex)
{
    std::vector<Eigen::Vector2i> face_sides_of_one_ring_vertices;
    vertex_one_ring_vertices(F,G,v2fs,v,one_ring_vertices,is_boundary_vertex,face_sides_of_one_ring_vertices);
}

void vertex_one_ring_vertices(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXi & v2fs,
    const int & v,
    Eigen::VectorXi & one_ring_vertices)
{
    bool is_boundary_vertex;
    std::vector<Eigen::Vector2i> face_sides_of_one_ring_vertices;
    vertex_one_ring_vertices(F,G,v2fs,v,one_ring_vertices,is_boundary_vertex,face_sides_of_one_ring_vertices);
}
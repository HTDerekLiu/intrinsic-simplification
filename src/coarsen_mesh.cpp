#include "coarsen_mesh.h"

static void initialize_transport_cost(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXi & G,
    const Eigen::MatrixXd & l,
    const Eigen::MatrixXi & v2fs,
    const Eigen::VectorXd & K,
    const double & w, // mixture weight (0~1) between curvature and area. 0 for pure curvature, 1 for pure area
    Eigen::MatrixXd & T)
{
    using namespace std;
    int nV = v2fs.rows();
    T.resize(nV, 9);
    T.setZero();

    // initialize cost 
    Eigen::VectorXd VA; // vertex area
    vertex_areas(F, l, VA);
    VA = VA / VA.sum(); // normalize to unit area
    double Kv; 
    for (int v=0; v<nV; v++)
    {
        // add curvature to the first two channels
        Kv = K(v);

        if (Kv >= 0)
            T(v, 0*3) = Kv * (1.-w);
        else
            T(v, 1*3) = -Kv * (1.-w);

        // add area to the third channel
        T(v, 2*3) = VA(v) * w;
    }
}

void coarsen_mesh(
    const int & total_removal,
    const double & mixture_weight, // 0 is pure curvature, 1 is pure area
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,    
    Eigen::MatrixXd & BC,
    std::vector<std::vector<int>> & F2V)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    int nV = v2fs.rows();
    int nF = F.rows();

    // compute Gaussian curvature and keep track of it (save some computation time)
    Eigen::VectorXd K(nV);
    for (int v=0; v<nV; v++)
        K(v) = gaussian_curvature_at_vertex(G,l,v2fs,v);

    // initialize transport cost
    Eigen::MatrixXd T;
    initialize_transport_cost(F, G,l,v2fs,K,mixture_weight,T);

    // initialize barycentric information
    if (BC.rows() < nV || BC.cols() < 3) { // if BC hasn't been initialized
        BC.resize(nV,3); 
        BC.setConstant(DOUBLE_INF);
    } else { // if BC has been initialized (usually for visualization)
        BC.block(0, 0, nV, 3).setConstant(DOUBLE_INF);
    }
    if (F2V.size() != nF)
        F2V.resize(nF);

    // create a priority queue sorted with vertex costs
    // Note: 
    // cost = get<0>(v_min_heap.top());
    // time_stamp = get<1>(v_min_heap.top());
    // vertex_index = get<2>(v_min_heap.top());
    typedef tuple<double, int, int> qtype; 
    priority_queue<qtype, vector<qtype>, greater<qtype> > v_min_heap;
    for (int v=0; v<nV; v++)
    {
        double cost;
        if (is_interior_vertex(G,v2fs,v))
            cost = pre_flatten_interior_vertex_and_cost(F,G,l,A,v2fs,T,K,v);
        else if (is_ear_vertex(G,v2fs,v))
            cost = pre_flatten_ear_vertex_and_cost(F,G,l,A,v2fs,T,K,v);
        else // regular boundary vertex
            cost = pre_flatten_boundary_vertex_and_cost(F,G,l,A,v2fs,T,K,v);
        qtype element(cost, 0, v);
        v_min_heap.push(element);
    }

    // create a time stamp vector for checking the time when v was inserted. This is useful to avoid updating min heap
    VectorXi v_time_stamp(nV);
    v_time_stamp.setZero();

    int num_removed = 0;
    double v_cost;
    while (v_min_heap.size() > 0 && num_removed<total_removal)
    {
        // get the top element
        qtype top_element = v_min_heap.top();
        v_min_heap.pop();    

        // STOP: if all vertices cannot be removed
        v_cost = get<0>(top_element);
        if (isinf(v_cost))
        {
            cout << "removed " << num_removed << "/" << total_removal << " vertices and it cannot be coarsened further" << endl;
            return;
        } 

        // get vertex index
        int v = get<2>(top_element);

        // SKIP: if time stamp is not the latest one
        int time_stamp = get<1>(top_element);
        if (time_stamp != v_time_stamp(v)) // if time stamp not match -> this element is an older element
            continue;

        int v_type;
        if (is_interior_vertex(G,v2fs,v))
            v_type = 0;
        else if (is_ear_vertex(G,v2fs,v))
            v_type = 1;
        else // regular boundary vertex
            v_type = 2;
        
        // extract one-ring for adding them into the queue later
        VectorXi v_onering;
        bool is_boundary_vertex;
        vertex_one_ring_vertices(F,G,v2fs,v,v_onering,is_boundary_vertex);

        // SKIP: if this vertex has been removed
        if (v_onering.size() == 0)
            continue;

        // flatten and compute the cost
        if (v_type == 0) // interior
            v_cost = flatten_interior_vertex_and_cost(F,G,v2fs,F2V,v,l,A,BC,T,K);
        else if (v_type == 1) // ear
            v_cost = flatten_ear_vertex_and_cost(F,G,v2fs,F2V,v,l,A,BC,T,K);
        else if (v_type == 2)  //regular boundary
            v_cost = flatten_boundary_vertex_and_cost(F,G,v2fs,F2V,v,l,A,BC,T,K);
        
        // SKIP: if flattening failed, skip this vertex. (intrinsic data won't be changed in place either)
        if (isinf(v_cost))
        {
            cout << "cetm failed at vertex " << v << ", skip it" << endl;
            continue;
        }

        // if flattening was successful, then go ahead and remove it
        if (v_type == 0) // interior
        {
            always_flip_to_degree_three(v,F,G,l,A,v2fs,BC,F2V);
            remove_degree_three_vertex(v,F,G,l,A,v2fs,BC,F2V,T);
        }
        else if (v_type == 1) // ear
        {
            remove_ear_vertex(v,F,G,l,A,v2fs,BC,F2V,T);
        }
        else if (v_type == 2) // regular boundary
        {
            always_flip_to_ear(v,F,G,l,A,v2fs,BC,F2V);
            remove_ear_vertex(v,F,G,l,A,v2fs,BC,F2V,T);
        }
        // after removal, count += 1
        num_removed += 1;
        if (num_removed % 5000 == 0)
            cout << "removed " << num_removed << " / " << total_removal << endl;

        // after removal, add one-ring vertices back to the queue
        for (int ii=0; ii<v_onering.size(); ii++)
        {
            int vv = v_onering(ii);
            if (v != vv){ // if this is not a self edge, then add vv back
                double cost;
                
                if (is_interior_vertex(G,v2fs,vv))
                    cost = pre_flatten_interior_vertex_and_cost(F,G,l,A,v2fs,T,K,vv);
                else if (is_ear_vertex(G,v2fs,vv))
                    cost = pre_flatten_ear_vertex_and_cost(F,G,l,A,v2fs,T,K,vv);
                else // regular boundary vertex
                    cost = pre_flatten_boundary_vertex_and_cost(F,G,l,A,v2fs,T,K,vv);

                qtype element(cost, num_removed, vv);
                v_min_heap.push(element);
                v_time_stamp(vv) = num_removed;
            }
        }
    }
}
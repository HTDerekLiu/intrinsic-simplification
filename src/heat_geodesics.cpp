#include "heat_geodesics.h"

bool heat_geodesics_precompute(
    const Eigen::MatrixXd & Vorl,
    const Eigen::MatrixXi & F,
    heat_geodesics_data & data)
{
    using namespace Eigen;
    using namespace std;

    // construct operators
    double t;
    SparseMatrix<double> L, M;
    VectorXd dblA;
    if (Vorl.rows() == F.rows()) // input is l, is intrinsic
    {
        igl::cotmatrix_intrinsic(Vorl,F,L);
        igl::massmatrix_intrinsic(Vorl,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
        igl::doublearea(Vorl,0,dblA);
        igl::grad_intrinsic(Vorl,F,data.Grad);

        // compute t
        const double avg_half_edge_length = Vorl.mean(); 
        t = avg_half_edge_length * avg_half_edge_length;
    }
    else // input is V, is extrinsic
    { 
        igl::cotmatrix(Vorl,F,L);
        igl::massmatrix(Vorl,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
        igl::doublearea(Vorl,F,dblA);
        igl::grad(Vorl,F,data.Grad);

        // compute t
        MatrixXd l;
        build_edge_lengths(Vorl,F,l);
        const double avg_half_edge_length = l.mean(); 
        t = avg_half_edge_length * avg_half_edge_length;
    }

    // div
    assert(F.cols() == 3 && "Only triangles are supported");
    // number of gradient components
    data.ng = data.Grad.rows() / F.rows();
    assert(data.ng == 3 || data.ng == 2);
    data.Div = -0.25*data.Grad.transpose()*dblA.colwise().replicate(data.ng).asDiagonal();

    Eigen::SparseMatrix<double> Q = M - t*L;
    Eigen::MatrixXi O;
    igl::boundary_facets(F,O);
    igl::unique(O,data.b);
    {
        Eigen::SparseMatrix<double> _;
        if(!igl::min_quad_with_fixed_precompute(
        Q,Eigen::VectorXi(),_,true,data.Neumann))
            return false;

        // Only need if there's a boundary
        if(data.b.size()>0)
        {
            if(!igl::min_quad_with_fixed_precompute(Q,data.b,_,true,data.Dirichlet))
                return false;
        }

        const RowVectorXd M_diag_tr = M.diagonal().transpose();
        const Eigen::SparseMatrix<double> Aeq = M_diag_tr.sparseView();
        L *= -0.5;
        if(!igl::min_quad_with_fixed_precompute(
        L,Eigen::VectorXi(),Aeq,true,data.Poisson))
            return false;
    }
    return true;
}

void heat_geodesics_solve(
    const heat_geodesics_data & data,
    const Eigen::VectorXi & gamma,
    Eigen::VectorXd & D)
{
    using namespace Eigen;

    // number of mesh vertices
    const int n = data.Grad.cols();
    // Set up delta at gamma
    VectorXd u0 = VectorXd::Zero(n,1);
    for(int g = 0;g<gamma.size();g++)
    {
        u0(gamma(g)) = 1;
    }
    // Neumann solution
    VectorXd u;
    igl::min_quad_with_fixed_solve(data.Neumann,u0,VectorXd(),VectorXd(),u);
    if(data.b.size()>0)
    {
        // Average Dirichelt and Neumann solutions
        VectorXd uD;
        igl::min_quad_with_fixed_solve(
        data.Dirichlet,u0,VectorXd::Zero(data.b.size()).eval(),VectorXd(),uD);
        u += uD;
        u *= 0.5;
    }
    VectorXd grad_u = data.Grad*u;
    const int m = data.Grad.rows()/data.ng;
    for(int i = 0;i<m;i++)
    {
        // It is very important to use a stable norm calculation here. If the
        // triangle is far from a source, then the floating point values in the
        // gradient can be _very_ small (e.g., 1e-300). The standard/naive norm
        // calculation will suffer from underflow. Dividing by the max value is more
        // stable. (Eigen implements this as stableNorm or blueNorm).
        double norm = 0;
        double ma = 0;
        for(int d = 0;d<data.ng;d++) {ma = std::max(ma,std::fabs(grad_u(d*m+i)));}
        for(int d = 0;d<data.ng;d++)
        {
        const double gui = grad_u(d*m+i) / ma;
        norm += gui*gui;
        }
        norm = ma*sqrt(norm);
        // These are probably over kill; ma==0 should be enough
        if(ma == 0 || norm == 0 || norm!=norm)
        {
        for(int d = 0;d<data.ng;d++) { grad_u(d*m+i) = 0; }
        }else
        {
        for(int d = 0;d<data.ng;d++) { grad_u(d*m+i) /= norm; }
        }
    }
    const VectorXd div_X = -data.Div*grad_u;
    const VectorXd Beq = (VectorXd(1,1)<<0).finished();
    igl::min_quad_with_fixed_solve(data.Poisson,(-div_X).eval(),VectorXd(),Beq,D);
    VectorXd Dgamma;
    igl::slice(D,gamma,Dgamma);
    D.array() -= Dgamma.mean();
    if(D.mean() < 0)
    {
        D = -D;
    }
}

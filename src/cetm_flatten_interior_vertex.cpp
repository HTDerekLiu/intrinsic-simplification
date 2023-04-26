#include "cetm_flatten_interior_vertex.h"

static void opposite_angles_from_lengths(
    const Eigen::VectorXd & lopps,
    const Eigen::VectorXd & lbs,
    const Eigen::VectorXd & lcs,
    Eigen::VectorXd & angles)
{
    using namespace std;
    int nF = lopps.size();
    angles.resize(nF);

    double lopp, lb, lc, lsum, scale, d;
    for (int ii=0; ii<nF; ii++)
    {
        lopp = lopps(ii);
        lb = lbs(ii);
        lc = lcs(ii);
        lsum = lopp + lb + lc;
        scale = 1.0 / lsum;
        lopp = lopp * scale;
        lb = lb * scale;
        lc = lc * scale;
        d = (lb*lb + lc*lc - lopp*lopp) / (2*lb*lc);
        angles(ii) = acos(d);
    }
}

static double gaussian_curvature_from_lengths(
    const Eigen::MatrixXd & l012)
{
    Eigen::VectorXd angles;
    opposite_angles_from_lengths(l012.col(1), l012.col(2), l012.col(0), angles);
    double angle_sum = angles.sum();

    // note that since the input vertex is an interior vertex so we only need to compute 2pi - anglesum
    return 2 * M_PI - angle_sum;
}

static double delta_u(
    const Eigen::MatrixXd & l012)
{
	Eigen::VectorXd l0 = l012.col(0);
	Eigen::VectorXd l1 = l012.col(1);
	Eigen::VectorXd l2 = l012.col(2);
    int nF = l012.rows();
    
    // gradient computation
    Eigen::VectorXd angles;
    opposite_angles_from_lengths(l1, l2, l0, angles);
    double gradient = 2*M_PI - angles.sum();

    // hessian computation
    Eigen::VectorXd alphas, cot_alphas;
    opposite_angles_from_lengths(l0, l1, l2, alphas);
    cot_alphas.resize(nF);
    for (int ii=0; ii<nF; ii++)
    {
        if (alphas(ii)<=0.0 || alphas(ii)>=M_PI)
            cot_alphas(ii) = 0.0;
        else
            cot_alphas(ii) = 1.0 / tan(alphas(ii));
    }
    Eigen::VectorXd betas, cot_betas;
    opposite_angles_from_lengths(l2, l0, l1, betas);
    cot_betas.resize(nF);
    for (int ii=0; ii<nF; ii++)
    {
        if (betas(ii)<=0.0 || betas(ii)>=M_PI)
            cot_betas(ii) = 0.0;
        else
            cot_betas(ii) = 1.0 / tan(betas(ii));
    }
    double hessian = (cot_betas.sum() + cot_alphas.sum()) / 4.0;

    return gradient * 0.5 / hessian;
}

static void scale_edge_lengths(
	const Eigen::MatrixXd & l012,
	const Eigen::MatrixXd & s012,
	const double & u, // conformal scaling factor
	Eigen::MatrixXd & l012_scaled)
{
	using namespace Eigen;
	using namespace std;
	MatrixXd u012 = s012.array() * u;
	MatrixXd scale012 = exp(u012.array() / 2.0);
	l012_scaled.resize(l012.rows(), l012.cols());
	l012_scaled = l012.array() * scale012.array();
}

void cetm_flatten_interior_vertex(
    const Eigen::MatrixXd & l012,
    const Eigen::MatrixXd & s012,
    double & u)
{
    using namespace std;
    using namespace Eigen;

    u = 0.0; // initial u
    double convergence_threshold = 1e-7;

    double Kv = gaussian_curvature_from_lengths(l012);
    int max_iter = 50;
    int decrease_iter = 10;
    double rho = 0.5; // decreasing factor for backtracking 
    double scale = 1.0; // actual scaling factor = exp(u/2)
    MatrixXd l012_scaled = l012; // initialized scaled edge lengths
    for (int iter=0; iter<max_iter; iter++)
    {
        if (abs(Kv) < convergence_threshold) break;

        // compute desent step using newton's method
        double du = delta_u(l012_scaled);

        // compute new curvature
        double step_size = 1.0;
		scale_edge_lengths(l012, s012, u-step_size*du, l012_scaled);
        Kv = gaussian_curvature_from_lengths(l012_scaled);

        // if NAN, backtracing
        for (int jj=0; jj<decrease_iter; jj++)
        {
            if (isnan(Kv))
            {
				// if NAN, scale down step size and compute new edge lengths
                step_size *= rho;
                scale_edge_lengths(l012, s012, u-step_size*du, l012_scaled);
            }
            else
            {
				// if not NAN, commit the change to the scaling factor u
                u = u - step_size * du;
                break;
            }
			// recompute Gaussian curvature from the newly scaled edge lengths
            Kv = gaussian_curvature_from_lengths(l012_scaled);
        }

        // after backtracing if it is still NAN
        if (isnan(Kv))
        {
            // cout << "CETM violates triangle inequalities" << endl;
            u = std::numeric_limits<double>::quiet_NaN();
            return;
        }
    }

    // if not converged
    if (abs(Kv) >= convergence_threshold)
    {
        // cout << "CETM did not converge" << endl;
        u = std::numeric_limits<double>::quiet_NaN();
        return;
    }
}

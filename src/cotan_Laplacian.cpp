#include "cotan_Laplacian.h"

void cotan_Laplacian(
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & l,
    Eigen::SparseMatrix<double> & L)
{   
    using namespace std;
    using namespace Eigen;

    int nV = F.maxCoeff() + 1;
    int nF = F.rows();

    vector<Triplet<double>> IJV;
    IJV.reserve(nF * 3 * 4);

    for (int f=0; f<F.rows(); f++){
        for (int s=0; s<3; s++){
            
            int i = F(f,s);
            int j = F(f,(s+1)%3);

            Vector2i fs; 
            fs << f, s;
            double opp_theta = opposite_corner_angle(l, fs);
            double opp_cotan =  1. / tan(opp_theta);
            double cotan_weight = 0.5 * opp_cotan;

            IJV.emplace_back(i, j, -cotan_weight);
            IJV.emplace_back(j, i, -cotan_weight);
            IJV.emplace_back(i, i, cotan_weight);
            IJV.emplace_back(j, j, cotan_weight);
        }
    }
    L.resize(nV,nV);
    L.setFromTriplets(IJV.begin(), IJV.end());
}
#include "get_prolongation.h"

void get_prolongation(
	const Eigen::MatrixXi & F, 
    const Eigen::MatrixXd & BC, 
    const std::vector<std::vector<int>> & F2V,
    const std::map<int, int> & IMV,
    Eigen::SparseMatrix<double> & P)
{
    using namespace Eigen;
    using namespace std;

    int nVi = BC.rows();
    int nV = F.maxCoeff() + 1;

    vector<Triplet<double>> IJV;
    IJV.reserve(3 * nV); // should be sufficient

    // add data for vertices were kept
    for (int vi=0; vi<nVi; vi++)
        if (isinf(BC(vi,0))) // meaning this vertex was kept
            IJV.push_back(Triplet<double>(vi, IMV.at(vi), 1.0));

    // add data for removed vertices
    int nF = F2V.size();
    for (int fIdx=0; fIdx<nF; fIdx++)
    {
        if (F2V[fIdx].size() > 0)
        {
            // gather vertex positions for this face
            int v0, v1, v2;
            v0 = F(fIdx,0);
            v1 = F(fIdx,1);
            v2 = F(fIdx,2);

            for (auto it = std::begin(F2V[fIdx]); it != std::end(F2V[fIdx]); ++it)
            {
                int idx = *it;
                IJV.push_back(Triplet<double>(idx, v0, BC(idx,0)));
                IJV.push_back(Triplet<double>(idx, v1, BC(idx,1)));
                IJV.push_back(Triplet<double>(idx, v2, BC(idx,2)));
            }
        }
    }
    P.resize(nVi, nV);
    P.setFromTriplets(IJV.begin(), IJV.end());
}
#include "get_barycentric_points.h"

void get_barycentric_points(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & BC,
    const std::vector<std::vector<int>> & F2V,
    Eigen::MatrixXd & P)
{
    using namespace Eigen;
    using namespace std;
    
    // BC may contain redundant barycentrics. Thus we keep track of "valid_indices" and slice the valid one to "P" later
    MatrixXd Ptmp(BC.rows(), 3);
	vector<int> valid_indices; 
	
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

			Vector3d p0, p1, p2;
			p0 << V.row(v0).transpose();
			p1 << V.row(v1).transpose();
			p2 << V.row(v2).transpose();

			for (auto it = std::begin(F2V[fIdx]); it != std::end(F2V[fIdx]); ++it)
			{
				int idx = *it;
				valid_indices.push_back(idx);
				
				double b0, b1, b2;
				b0 = BC(idx,0);
				b1 = BC(idx,1);
				b2 = BC(idx,2);

				Ptmp.row(idx) = (b0*p0 + b1*p1 + b2*p2).transpose();
			}
			// std::sort(valid_indices.begin(), valid_indices.end());
		}
	}
	std::sort(valid_indices.begin(), valid_indices.end());

    // slice Ptmp back to P
	int num_valid_indices = valid_indices.size();
	P.resize(num_valid_indices,3);
    P.setZero();
	int pIdx = 0;
	for (int ii=0; ii<num_valid_indices; ii++)
	{
		int idx = valid_indices[ii];
		P.row(pIdx) = Ptmp.row(idx);
		pIdx += 1;
	}
}

void get_barycentric_points_redundant(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & BC,
    const std::vector<std::vector<int>> & F2V,
    Eigen::MatrixXd & P)
{
    using namespace Eigen;
    using namespace std;
    
    // BC may contain redundant barycentrics. Thus we keep track of "valid_indices" and slice the valid one to "P" later
    P.resize(BC.rows(), 3);
	P.setZero();
	
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

			Vector3d p0, p1, p2;
			p0 << V.row(v0).transpose();
			p1 << V.row(v1).transpose();
			p2 << V.row(v2).transpose();

			for (auto it = std::begin(F2V[fIdx]); it != std::end(F2V[fIdx]); ++it)
			{
				int idx = *it;
				double b0, b1, b2;
				b0 = BC(idx,0);
				b1 = BC(idx,1);
				b2 = BC(idx,2);
				P.row(idx) = (b0*p0 + b1*p1 + b2*p2).transpose();
			}
		}
	}
}
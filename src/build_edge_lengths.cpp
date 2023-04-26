#include "build_edge_lengths.h"

void build_edge_lengths(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & l)
{
    int nF = F.rows();
	l.resize(nF,3);
	for (int f=0; f<nF; f++)
	{
		for (int s=0; s<3; s++)
		{
			int vi = F(f,s);
			int vj = F(f,(s+1)%3);
			l(f,s) = (V.row(vj) - V.row(vi)).norm();
		}
	}
}
#include "mg_precompute.h"

void mg_precompute(
    const Eigen::MatrixXi & Fi,
    const Eigen::MatrixXi & Gi,
    const Eigen::MatrixXd & li,
    const Eigen::MatrixXd & Ai,
    const Eigen::MatrixXi & v2fsi, 
	const double & ratio,
	const int & nV_coarsest,
	const double & mixture_weight,
  	std::vector<mg_data> & mg)
{
	using namespace std;
	using namespace Eigen;

	// compute number of multigrid levels
	int nLvs;
	{
		nLvs = 1;
		int nV = Fi.maxCoeff() + 1;
		while (nV > nV_coarsest)
		{
			nV *= ratio;
            nLvs += 1;
		}
	}

    // assign the first mg level (TODO: do not save the information that is not needed)
	mg_data data_lv0;
    data_lv0.F = Fi;

	// create vertex index (first level is simply [0,nV-1])
	{
		VectorXi vIdx(Fi.maxCoeff() + 1);
		for (int v=0; v<vIdx.size(); v++)
			vIdx(v) = v;
		data_lv0.vIdx = vIdx;
	}

	mg.clear();
	mg.reserve(nLvs);
	mg.push_back(data_lv0);

	// initialize intrinsic information
	MatrixXi F = Fi;
	MatrixXi G = Gi;
	MatrixXi v2fs = v2fsi;
	MatrixXd l = li;
	MatrixXd A = Ai;
	for (int lv = 1; lv < nLvs; lv++)
	{
        int nV = mg[lv-1].F.maxCoeff()+1;
		int nV_tar = round((float)(nV) * ratio);

        // coarsen the mesh
        int total_removal = nV - nV_tar;
        MatrixXd BC;
        vector<vector<int>> F2V;
        coarsen_mesh(total_removal,mixture_weight, F,G,l,A,v2fs,BC,F2V);

		// remove unreferenced
		map<int, int> IMV, IMF;
		VectorXi vIdx;
		VectorXi fIdx;
		remove_unreferenced_intrinsic(F,G,l,A,v2fs,F2V,IMV,IMF,vIdx,fIdx);

        // build prolongation
		SparseMatrix<double> P;
		get_prolongation(F,BC,F2V,IMV,P);
		cout << "lv: " << lv << ", Vc: " << vIdx.size() << endl;

		// create vIdx_c
		VectorXi vIdx_c(vIdx.size());
		VectorXi vIdx_f = mg[lv-1].vIdx;
		for (int ii=0; ii<vIdx.size(); ii++)
			vIdx_c(ii) = vIdx_f(vIdx(ii));

		mg_data data;
		data.F = F;
		data.vIdx = vIdx_c;

		data.P = P;
		data.PT = P.transpose();
		data.P_full = P;
		mg.push_back(data);
	}

	// print multigrid info
	cout << "============\n";
	cout << "Multigrid Info\n";
	cout << "============\n";
	cout << "numLv: " << mg.size() << endl;
	cout << "|V_coarsest|: " << mg[mg.size()-1].F.maxCoeff()+1 << endl;
}

void mg_precompute(
    const Eigen::MatrixXi & Fi,
    const Eigen::MatrixXi & Gi,
    const Eigen::MatrixXd & li,
    const Eigen::MatrixXd & Ai,
    const Eigen::MatrixXi & v2fsi, 
	const double & ratio,
	const int & nV_coarsest,
  	std::vector<mg_data> & mg)
{
	double mixture_weight = 0.5;
	mg_precompute(Fi,Gi,li,Ai,v2fsi,ratio,nV_coarsest,mixture_weight,mg);
}
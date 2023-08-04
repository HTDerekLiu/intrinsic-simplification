#ifndef MG_DATA
#define MG_DATA

#include <vector>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

struct mg_data
{
    // intrinsic data
	Eigen::MatrixXi F; // faces
	Eigen::VectorXi vIdx; // index of the remaining vertices (from the orignal V)

    // multigrid data
	Eigen::SparseMatrix<double> P_full; // full prolongation
	Eigen::SparseMatrix<double> A; // LHS for only unknown parts
	Eigen::VectorXd A_diag; // diagonal entries of A
	Eigen::SparseMatrix<double> P; // prolongation for only unknown parts 
	Eigen::SparseMatrix<double> PT; // prolongation transpose for only unknown parts 

	// Gauss Seidel relaxation precomputation
	std::vector<std::vector<int>> S;
	Eigen::VectorXi SV;
	Eigen::VectorXi SVI;
	Eigen::VectorXi SC;
	Eigen::VectorXi SCS;

	void reset()
	{
		F = Eigen::MatrixXi();
		vIdx = Eigen::VectorXi();

		P_full = Eigen::SparseMatrix<double>();
		A = Eigen::SparseMatrix<double>();
		A_diag = Eigen::VectorXd();
		// Ar = Eigen::SparseMatrix<double, Eigen::RowMajor>();
		P = Eigen::SparseMatrix<double>();
		PT = Eigen::SparseMatrix<double>();
        
		S.clear(); // reset the 2d array
		SV = Eigen::VectorXi();
		SVI = Eigen::VectorXi();
		SC = Eigen::VectorXi();
		SCS = Eigen::VectorXi();
	}
};

#endif
#include "remove_unreferenced_intrinsic.h"

void remove_unreferenced_intrinsic(
    Eigen::MatrixXi & F,
    std::map<int, int> & IMV,
    std::map<int, int> & IMF,
    Eigen::VectorXi & vIdx,
    Eigen::VectorXi & fIdx)
{    
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    // extract faces to keep
    vector<int> fIdx_vec;
    for (int f=0; f<F.rows(); f++)
        if (F(f,0) != GHOST_INDEX)
            fIdx_vec.push_back(f);

    // turn it into Eigen::VectorXi
    fIdx.resize(fIdx_vec.size());
    fIdx = Map<VectorXi, Unaligned>(fIdx_vec.data(), fIdx_vec.size());

    // assemble an index map IMF[old_fIdx] = new_fIdx
    IMF.clear();
    for(int ii = 0; ii < fIdx.size(); ii++)
        IMF[fIdx(ii)] = ii;

    // create a set of unique elements
    unordered_set<int> vIdx_set;
    for (int c=0; c<F.cols(); c++)
        for (int r=0; r<F.rows(); r++)
        {
            if (F(r,c) != GHOST_INDEX)
                vIdx_set.insert(F(r,c));
        }
    
    // sort the unique elements
    vector<int> vIdx_vec;
    vIdx_vec.assign( vIdx_set.begin(), vIdx_set.end() );
    sort( vIdx_vec.begin(), vIdx_vec.end() );

    // turn it into Eigen::VectorXi
    vIdx.resize(vIdx_vec.size());
    vIdx = Map<VectorXi, Unaligned>(vIdx_vec.data(), vIdx_vec.size());

    // assemble an index map IMV[old_vIdx] = new_vIdx
    IMV.clear();
    for(int ii = 0; ii < vIdx.size(); ii++)
        IMV[vIdx(ii)] = ii;

    // remove unreferenced F in place by overwriting it
    for (int c=0; c<F.cols(); c++)
        for (int r=0; r<fIdx.size(); r++)
        {
            int f = fIdx(r);
            F(r,c) = IMV.at(F(f,c));
        }
    F.conservativeResize(fIdx.size(),F.cols());
}

void remove_unreferenced_intrinsic(
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & G,
    Eigen::MatrixXd & l,
    Eigen::MatrixXd & A,
    Eigen::MatrixXi & v2fs,  
    std::vector<std::vector<int>> & F2V,  
    std::map<int, int> & IMV,
    std::map<int, int> & IMF,
    Eigen::VectorXi & vIdx,
    Eigen::VectorXi & fIdx)
{
    using namespace std;
    using namespace Eigen;
    using namespace global_variables;

    remove_unreferenced_intrinsic(F,IMV,IMF,vIdx,fIdx);
    int nF = fIdx.size();
    int nV = vIdx.size();

    // update G, l, A
    for (int c=0; c<3; c++)
        for (int r=0; r<nF; r++)
        {
            int f = fIdx(r);

            // update G
            if (G(f,c*2) != GHOST_INDEX)
            {
                G(r,c*2) = IMF.at(G(f,c*2));
                G(r,c*2+1) = G(f,c*2+1);
            }
            else
            {
               G(r,c*2) = GHOST_INDEX;
               G(r,c*2+1) = GHOST_INDEX;
            }

            // update l
            l(r,c) = l(f,c);

            // update A
            A(r,c) = A(f,c);
        }
    G.conservativeResize(nF,6);
    l.conservativeResize(nF,3);
    A.conservativeResize(nF,3);

    // update v2fs
    for (int ii=0; ii<nV; ii++)
    {   
        int v = vIdx(ii);
        v2fs(ii,0) = IMF.at(v2fs(v,0));
        v2fs(ii,1) = v2fs(v,1);
    }
    v2fs.conservativeResize(nV,2);

    // remove lists from F2V by overwriting and then pop back
    for (int ii=0; ii<nF; ii++)
    {
        int f = fIdx(ii);
        F2V[ii] = F2V[f];
    }
    F2V.resize(nF);
}
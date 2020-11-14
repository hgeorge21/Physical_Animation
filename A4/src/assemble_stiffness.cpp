#include <assemble_stiffness.h>
#include <iostream>
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) { 
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
       
    for (int i = 0; i < F.rows(); i++) {
        Eigen::Matrix99d H;
        Eigen::RowVector3i element = F.row(i);

        Eigen::Matrix<double, 1, 9> tmp_row;
        tmp_row = dX.row(i);
        Eigen::Map<const Eigen::Matrix3d> dXi(tmp_row.data());
        
        d2V_membrane_corotational_dq2(H, q, dXi, V, element, a0(i), mu, lambda);

        for (int l = 0; l < 3; l++) {
            for (int k = 0; k < 3; k++) {
                triplets.push_back(T(3 * F(i, l), 3 * F(i, k),     H(3 * l, 3 * k)));
                triplets.push_back(T(3 * F(i, l), 3 * F(i, k)+1,   H(3 * l, 3 * k+1)));
                triplets.push_back(T(3 * F(i, l), 3 * F(i, k)+2,   H(3 * l, 3 * k+2)));
                triplets.push_back(T(3 * F(i, l)+1, 3 * F(i, k),   H(3 * l+1, 3 * k)));
                triplets.push_back(T(3 * F(i, l)+1, 3 * F(i, k)+1, H(3 * l+1, 3 * k+1)));
                triplets.push_back(T(3 * F(i, l)+1, 3 * F(i, k)+2, H(3 * l+1, 3 * k+2)));
                triplets.push_back(T(3 * F(i, l)+2, 3 * F(i, k),   H(3 * l+2, 3 * k)));
                triplets.push_back(T(3 * F(i, l)+2, 3 * F(i, k)+1, H(3 * l+2, 3 * k+1)));
                triplets.push_back(T(3 * F(i, l)+2, 3 * F(i, k)+2, H(3 * l+2, 3 * k+2)));
            }
        }
    }
    
    K.resize(q.rows(), q.rows());
    K.setZero();
    K.setFromTriplets(triplets.begin(), triplets.end());
    K = -1. * K;
};

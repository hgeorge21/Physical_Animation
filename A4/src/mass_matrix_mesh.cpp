#include <mass_matrix_mesh.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {
    double h = 1.;
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    
    Eigen::Matrix3d I1 = Eigen::Matrix3d::Identity() / 12;
    Eigen::Matrix3d I2 = Eigen::Matrix3d::Identity() / 24;
    Eigen::MatrixXd mass;
    mass.resize(9, 9);
    mass << I1, I2, I2,
        I2, I1, I2,
        I2, I2, I1;
    mass = h * density * mass;

    for (int i = 0; i < F.rows(); i++) {
        Eigen::MatrixXd tmp;
        tmp = 2. * areas(i) * mass;

        for (int l = 0; l < 3; l++) {
            for (int k = 0; k < 3; k++) {
                triplets.push_back(T(3 * F(i, l), 3 * F(i, k),   tmp(3*l, 3*k)));
                triplets.push_back(T(3 * F(i, l), 3 * F(i, k)+1, tmp(3*l, 3*k+1)));
                triplets.push_back(T(3 * F(i, l), 3 * F(i, k)+2, tmp(3*l, 3*k+2)));
                triplets.push_back(T(3 * F(i, l)+1, 3 * F(i, k),   tmp(3*l+1, 3*k)));
                triplets.push_back(T(3 * F(i, l)+1, 3 * F(i, k)+1, tmp(3*l+1, 3*k+1)));
                triplets.push_back(T(3 * F(i, l)+1, 3 * F(i, k)+2, tmp(3*l+1, 3*k+2)));
                triplets.push_back(T(3 * F(i, l)+2, 3 * F(i, k),   tmp(3*l+2, 3*k)));
                triplets.push_back(T(3 * F(i, l)+2, 3 * F(i, k)+1, tmp(3*l+2, 3*k+1)));
                triplets.push_back(T(3 * F(i, l)+2, 3 * F(i, k)+2, tmp(3*l+2, 3*k+2)));
            }
        }
    }
   
    M.resize(3 * V.rows(), 3 * V.rows());
    M.setFromTriplets(triplets.begin(), triplets.end());
}
 
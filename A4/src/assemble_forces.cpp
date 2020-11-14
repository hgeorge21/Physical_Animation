#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) {
    f.resize(q.rows());
    f.setZero();
    for (int i = 0; i < F.rows(); i++) {
        Eigen::Vector9d force;
        Eigen::Matrix<double, 1, 9> tmp_row;
        tmp_row = dX.row(i);
        Eigen::Map<const Eigen::Matrix3d> dXi(tmp_row.data());
        
        dV_membrane_corotational_dq(force, q, dXi, V, F.row(i), a0(i), mu, lambda);
                    
        for (int j = 0; j < 3; j++) {
            f(3 * F(i, j)) += force(3 * j);
            f(3 * F(i, j) + 1) += force(3 * j + 1);
            f(3 * F(i, j) + 2) += force(3 * j + 2);
        }
    }     
    f = -1. * f;
};

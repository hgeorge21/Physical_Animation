#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D) { 
    f.resize(3 * q.rows());
    f.setZero();
    for (int i = 0; i < T.rows(); i++) {
        Eigen::Vector12d fi;
        Eigen::RowVector4i element = T.row(i);

        double volume = v0(i);
        dV_linear_tetrahedron_dq(fi, q, V, element, volume, C, D);

        for (int j = 0; j < 4; j++) {
            f(3 * T(i, j))   -= fi(3 * j);
            f(3 * T(i, j)+1) -= fi(3 * j+1);
            f(3 * T(i, j)+2) -= fi(3 * j+2);
        }
    }
};
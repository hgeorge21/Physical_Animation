#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    Eigen::Vector3d X0 = V.row(element(0));
    Eigen::Vector3d X1 = V.row(element(1));
    Eigen::Vector3d X2 = V.row(element(2));
    Eigen::Vector3d X3 = V.row(element(3));

    Eigen::MatrixXd T(3, 3);
    T << X1 - X0, X2 - X0, X3 - X0;

    Eigen::MatrixXd D(4, 3);
    D.block<1, 3>(0, 0) = -1 * Eigen::RowVector3d::Ones();
    D.block<3, 3>(1, 0) = Eigen::Matrix3d::Identity();

    dphi = D * T.inverse();
}
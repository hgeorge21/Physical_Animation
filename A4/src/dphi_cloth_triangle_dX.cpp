#include <dphi_cloth_triangle_dX.h>
#include <iostream>
//compute 3x3 deformation gradient 
void dphi_cloth_triangle_dX(Eigen::Matrix3d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    Eigen::MatrixXd T(3, 2);
    T.block<3, 1>(0, 0) = V.row(element(1)) - V.row(element(0));
    T.block<3, 1>(0, 1) = V.row(element(2)) - V.row(element(0));

    Eigen::MatrixXd temp = (T.transpose() * T).inverse() * T.transpose();

    Eigen::Matrix32d tmp;
    tmp << -1. * Eigen::RowVector2d::Ones(), 
            Eigen::Matrix2d::Identity();
    dphi = tmp * temp;
}
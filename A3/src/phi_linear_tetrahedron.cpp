#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {
    Eigen::Vector3d x0 = V.row(element(0));
    Eigen::Vector3d x1 = V.row(element(1));
    Eigen::Vector3d x2 = V.row(element(2));
    Eigen::Vector3d x3 = V.row(element(3));
    
    Eigen::MatrixXd T(3, 3);
    T << x1-x0, x2-x0, x3-x0;

    Eigen::Vector3d tmp = T.inverse() * (x - x0);
    double phi0 = 1 - tmp.sum();

    phi[0] = phi0;
    phi.segment<3>(1) = tmp;
}
#include <rigid_body_jacobian.h>

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> X) {
    Eigen::Vector3d x = X - p;
    Eigen::Matrix3d Xbar;
    Xbar << 0, -x(2), x(1),
            x(2), 0, -x(0),
            -x(1), x(0), 0;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d O = Eigen::Matrix3d::Zero();

    Eigen::Matrix36d tmp1;
    tmp1 << Xbar.transpose(), I;

    Eigen::Matrix66d tmp2;
    tmp2 << R.transpose(), O,
           O, I;
    
    J = R * tmp1 * tmp2;
}


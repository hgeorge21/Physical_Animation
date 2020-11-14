#include <V_membrane_corotational.h>
#include <dphi_cloth_triangle_dX.h>

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    Eigen::Vector3d x0 = q.segment<3>(3 * element(0));
    Eigen::Vector3d x1 = q.segment<3>(3 * element(1));
    Eigen::Vector3d x2 = q.segment<3>(3 * element(2));
    Eigen::Vector3d n = (x1 - x0).cross(x2 - x0); // .normalized();

    Eigen::Vector3d X0 = V.row(element(0));
    Eigen::Vector3d X1 = V.row(element(1));
    Eigen::Vector3d X2 = V.row(element(2));
    Eigen::Vector3d N = (X1 - X0).cross(X2 - X0); // .normalized();

    Eigen::Matrix34d tmp1;
    tmp1 << x0, x1, x2, n;

    Eigen::Matrix43d tmp2;
    tmp2 << dX, N.transpose();

    Eigen::Matrix3d F = tmp1 * tmp2;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d S = svd.singularValues();

    energy = 0;
    for (int i = 0; i < 3; i++)
        energy += mu * pow(S(i) - 1, 2);    
    energy += 0.5 * lambda * pow(S.sum() - 3, 2);
}

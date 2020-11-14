#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>
void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D) {
    double J = F.determinant();
    double I = (F.transpose() * F).trace();

    psi = C * (I / pow(J, 2. / 3) - 3) + D * pow(J - 1, 2);
}
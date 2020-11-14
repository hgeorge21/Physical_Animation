#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>
#include <iostream>

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        Eigen::MatrixXd tmpX(3, 4);
        tmpX.block<3, 1>(0, 0) = q.segment<3>(3 * element(0));
        tmpX.block<3, 1>(0, 1) = q.segment<3>(3 * element(1));
        tmpX.block<3, 1>(0, 2) = q.segment<3>(3 * element(2));
        tmpX.block<3, 1>(0, 3) = q.segment<3>(3 * element(3));

        Eigen::Matrix43d dphi;
        dphi_linear_tetrahedron_dX(dphi, V, element, X);

        Eigen::Matrix3d F = tmpX * dphi;
        psi_neo_hookean(e, F, C, D);
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  
    
}
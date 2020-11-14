#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
       Eigen::MatrixXd tmpX(3, 4);
       tmpX.block<3, 1>(0, 0) = q.segment<3>(3 * element(0));
       tmpX.block<3, 1>(0, 1) = q.segment<3>(3 * element(1));
       tmpX.block<3, 1>(0, 2) = q.segment<3>(3 * element(2));
       tmpX.block<3, 1>(0, 3) = q.segment<3>(3 * element(3));

       Eigen::Matrix43d dphi;
       dphi_linear_tetrahedron_dX(dphi, V, element, X);
               
       Eigen::Matrix3d F = tmpX * dphi;
       Eigen::Vector9d dpsi;
       dpsi_neo_hookean_dF(dpsi, F, C, D);

       Eigen::MatrixXd B = Eigen::MatrixXd::Zero(9, 12);
       for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 4; j++) {
               B(i, 3 * j) = dphi(j, i);
               B(i+3, 3 * j+1) = dphi(j, i);
               B(i+6, 3 * j+2) = dphi(j, i);
           }
       }
       dV = B.transpose() * dpsi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
    
}
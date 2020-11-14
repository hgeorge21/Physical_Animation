#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>
#include <iostream>

void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

   auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
       Eigen::MatrixXd tmpX(3, 4);
       tmpX.block<3, 1>(0, 0) = q.segment<3>(3 * element(0));
       tmpX.block<3, 1>(0, 1) = q.segment<3>(3 * element(1));
       tmpX.block<3, 1>(0, 2) = q.segment<3>(3 * element(2));
       tmpX.block<3, 1>(0, 3) = q.segment<3>(3 * element(3));

       Eigen::Matrix43d dphi;
       dphi_linear_tetrahedron_dX(dphi, V, element, tmpX);
      
       Eigen::Matrix3d F = tmpX * dphi;
       Eigen::Matrix99d d2psi;
       d2psi_neo_hookean_dF2(d2psi, F, C, D);

       Eigen::MatrixXd B = Eigen::MatrixXd::Zero(9, 12);
       for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 4; j++) {
               B(i, 3 * j) = dphi(j, i);
               B(i + 3, 3 * j + 1) = dphi(j, i);
               B(i + 6, 3 * j + 2) = dphi(j, i);
           }
       }
       dV = B.transpose() * d2psi * B;
    };

    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);  

     // hard coded for tet, need to change size for hex
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 12; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();

}
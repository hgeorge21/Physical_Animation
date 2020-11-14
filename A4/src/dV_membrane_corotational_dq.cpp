#include <dV_membrane_corotational_dq.h>
#include <dsvd.h>
#include <iostream>

void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix 
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 

    //TODO: SVD Here
    Eigen::Vector3d x0 = q.segment<3>(3 * element(0));
    Eigen::Vector3d x1 = q.segment<3>(3 * element(1));
    Eigen::Vector3d x2 = q.segment<3>(3 * element(2));
    Eigen::Vector3d dx1 = x1 - x0;
    Eigen::Vector3d dx2 = x2 - x0;
    Eigen::Vector3d n = (dx1).cross(dx2).normalized();

    Eigen::Vector3d X0 = V.row(element(0));
    Eigen::Vector3d X1 = V.row(element(1));
    Eigen::Vector3d X2 = V.row(element(2));
    Eigen::Vector3d N = (X1 - X0).cross(X2 - X0).normalized();

    Eigen::Matrix34d tmp1;
    tmp1 << x0, x1, x2, n;

    Eigen::Matrix43d D;
    D << dX, N.transpose();

    Eigen::Matrix3d F = tmp1 * D;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    W = svd.matrixV();
    S = svd.singularValues();

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }
    
    //TODO: energy model gradient 
    double h = 1.0;
    Eigen::Matrix3d dS = Eigen::Matrix3d::Zero();
    
    for (int i = 0; i < 3; i++)
        dS(i, i) = 2 * mu * (S(i) - 1) + lambda * (S.sum() - 3);

    // compute dF/dq
    Eigen::MatrixXd NN = Eigen::MatrixXd::Zero(9, 3);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(9, 9);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            B(i, 3 * j) = D(j, i);
            B(i + 3, 3 * j + 1) = D(j, i);
            B(i + 6, 3 * j + 2) = D(j, i);
        }
        NN(i, 0) = N(i);
        NN(i + 3, 1) = N(i);
        NN(i + 6, 2) = N(i);
    }
        
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d O = Eigen::Matrix3d::Zero();
    Eigen::Matrix39d res1, res2;
    Eigen::Matrix3d tt;

    Eigen::Matrix39d tmp;
    tmp << -I, O, I; 
    tt << 0, -dx1(2), dx1(1),
        dx1(2), 0, -dx1(0),
        -dx1(1), dx1(0), 0;
    res1 = tt * tmp;

    tmp << -I, I, O;
    tt << 0, -dx2(2), dx2(1),
        dx2(2), 0, -dx2(0),
        -dx2(1), dx2(0), 0;
    res2 = tt * tmp;

    Eigen::Matrix39d dndq = 1. / (dx1).cross(dx2).norm() * (I - n * n.transpose()) * (res1 - res2);
    Eigen::MatrixXd dFdq = B + NN * dndq;

    Eigen::Matrix3d dV_ =  U * dS * W.transpose();
    dV.segment<3>(0) = dV_.row(0);
    dV.segment<3>(3) = dV_.row(1);
    dV.segment<3>(6) = dV_.row(2);
    dV = h * area * dFdq.transpose() * dV;
}

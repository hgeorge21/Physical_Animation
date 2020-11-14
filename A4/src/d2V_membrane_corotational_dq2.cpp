#include <d2V_membrane_corotational_dq2.h>
#include <iostream>


void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    

    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    Eigen::Matrix3d F; //deformation gradient
    
    double tol = 1e-5;
    
    //Compute SVD of F here
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

    F = tmp1 * D;
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    W = svd.matrixV();
    S = svd.singularValues();
   
    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {
        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
    }
    
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

    //TODO: compute H, the hessian of the corotational energy
    Eigen::Tensor3333d dU;
    Eigen::Tensor3333d dW;
    Eigen::Tensor333d dS;
    Eigen::MatrixXd d2psidF = Eigen::MatrixXd::Zero(9, 9);
    dsvd(dU, dS, dW, F);
    
    Eigen::Matrix3d dpsi = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d d2psi = lambda * Eigen::Matrix3d::Ones();
   
    for (int i = 0; i < 3; i++) {
        dpsi(i, i) = 2 * mu * (S(i) - 1) + lambda * (S.sum() - 3);
        d2psi(i, i) += 2 * mu;
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Eigen::Vector3d diag = d2psi * dS[i][j];
            Eigen::Matrix3d diagM = diag.asDiagonal();

            Eigen::Matrix3d tmp = Eigen::Matrix3d::Zero();
            tmp += dU[i][j] * dpsi * W.transpose();
            tmp += U * diagM * W.transpose();
            tmp += U * dpsi * dW[i][j].transpose();

            d2psidF.block<1, 3>(3 * i + j, 0) = tmp.row(0);
            d2psidF.block<1, 3>(3 * i + j, 3) = tmp.row(1);
            d2psidF.block<1, 3>(3 * i + j, 6) = tmp.row(2);
        }
    }

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

    Eigen::Matrix39d dndq = 1. / (dx1).cross(dx2).norm() * (I - n * n.transpose()) * (res1-res2); 
    Eigen::MatrixXd dFdq = B + NN * dndq;

    // compute H    
    double h = 1.0;
    H = dFdq.transpose() * d2psidF * dFdq;
    H = h * area * H;

    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    
    H = Evec * DiagEval * Evec.transpose();
    
}

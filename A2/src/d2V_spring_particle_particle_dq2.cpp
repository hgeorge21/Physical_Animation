#include <d2V_spring_particle_particle_dq2.h>
#include <iostream>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    Eigen::MatrixXd B(3, 6);
    B << -1. * Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Identity();

    Eigen::MatrixXd t0 = q1 - q0;
    double t1 = t0.norm();
    Eigen::MatrixXd T2 = B.transpose() * t0 * t0.transpose() * B;
    double t3 = stiffness * (t1 - l0);

    Eigen::MatrixXd tmp;
    H = -1.0 * ((stiffness / std::pow(t1, 2) - t3 / std::pow(t1, 3)) * T2 + t3 / t1 * B.transpose() * B).eval();
}
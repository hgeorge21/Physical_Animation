#include <dV_spring_particle_particle_dq.h>
#include <iostream>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    Eigen::Vector3d t0 = q1-q0;
    double t1 = t0.norm();
    
    f.segment<3>(0) = q0 - q1;
    f.segment<3>(3) = q1 - q0;
    f = stiffness * (t1 - l0) / t1 * f;
}
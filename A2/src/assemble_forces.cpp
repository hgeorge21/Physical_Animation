#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
    f = Eigen::VectorXd::Zero(q.size());

    for (int i = 0; i < E.rows(); i++) {
        Eigen::Vector3d v1, v2;
        v1 = V.row(E(i, 0));
        v2 = V.row(E(i, 1));

        auto li = l0(i);

        Eigen::Vector6d force;
        dV_spring_particle_particle_dq(force, v1, v2, li, k);

        f.segment<3>(3 * E(i, 0)) -= force.segment<3>(0);
        f.segment<3>(3 * E(i, 1)) += force.segment<3>(3);
    }
};

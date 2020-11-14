#include <vector>
#include <assemble_stiffness.h>
#include <iostream>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;

    for (int i = 0; i < E.rows(); i++) { 
        Eigen::Vector3d q0, q1;
        q0 = q.segment<3>(3 * E(i, 0));
        q1 = q.segment<3>(3 * E(i, 1));
        
        Eigen::Matrix66d Hi;
        d2V_spring_particle_particle_dq2(Hi, q0, q1, l0(i), k);

        for (int l = 0; l < 3; l++)
            for (int j = 0; j < 3; j++)
                triplets.push_back(T(3 * E(i, 0) + l, 3 * E(i, 0) + j, Hi(l, j)));
        for (int l = 0; l < 3; l++)
            for (int j = 0; j < 3; j++)
                triplets.push_back(T(3 * E(i, 0) + l, 3 * E(i, 1) + j, Hi(l, 3+j)));

        for (int l = 0; l < 3; l++)
            for (int j = 0; j < 3; j++)
                triplets.push_back(T(3 * E(i, 1) + l, 3 * E(i, 0) + j, Hi(3+j, l)));
        for (int l = 0; l < 3; l++)
            for (int j = 0; j < 3; j++)
                triplets.push_back(T(3 * E(i, 1) + l, 3 * E(i, 1) + j, Hi(3+l, 3+j)));
    }

    K.resize(q.size(), q.size());
    K.setZero();
    K.setFromTriplets(triplets.begin(), triplets.end());
};
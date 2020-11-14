#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {
    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> triplets;

    for (int l = 0; l < V_skin.rows(); l++) {
        int tet_idx = 0;
        double phi_norm = std::numeric_limits<double>::max();
        Eigen::Vector4d phi;

        for (int i = 0; i < T.rows(); i++) {
            Eigen::Vector3d x = V_skin.row(l);
            Eigen::Vector4d phi_;
            Eigen::RowVectorXi element = T.row(i);
            phi_linear_tetrahedron(phi_, V, element, x);
            
            if (phi_.norm() < phi_norm) {
                tet_idx = i;
                phi_norm = phi_.norm();
                phi = phi_;
            }
        }

        triplets.push_back(Triplet(l, T(tet_idx, 0), phi(0)));
        triplets.push_back(Triplet(l, T(tet_idx, 1), phi(1)));
        triplets.push_back(Triplet(l, T(tet_idx, 2), phi(2)));
        triplets.push_back(Triplet(l, T(tet_idx, 3), phi(3)));
    }

    N.resize(3 * V_skin.rows(), 3 * V.rows());
    N.setFromTriplets(triplets.begin(), triplets.end());
}
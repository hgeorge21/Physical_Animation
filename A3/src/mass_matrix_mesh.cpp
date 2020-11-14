#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {
    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> triplets;

    for (int i = 0; i < T.rows(); i++) {
        double volume = v0(i);
        Eigen::RowVector4i element = T.row(i);

        Eigen::Matrix1212d m;
        mass_matrix_linear_tetrahedron(m, qdot, element, density, volume);
    
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                triplets.push_back(Triplet(3 * T(i, j),   3 * T(i, k),   m(3*j,   3*k  )));
                triplets.push_back(Triplet(3 * T(i, j),   3 * T(i, k)+1, m(3*j,   3*k+1)));
                triplets.push_back(Triplet(3 * T(i, j),   3 * T(i, k)+2, m(3*j,   3*k+2)));
                triplets.push_back(Triplet(3 * T(i, j)+1, 3 * T(i, k),   m(3*j+1, 3*k  )));
                triplets.push_back(Triplet(3 * T(i, j)+1, 3 * T(i, k)+1, m(3*j+1, 3*k+1)));
                triplets.push_back(Triplet(3 * T(i, j)+1, 3 * T(i, k)+2, m(3*j+1, 3*k+2)));
                triplets.push_back(Triplet(3 * T(i, j)+2, 3 * T(i, k),   m(3*j+2, 3*k  )));
                triplets.push_back(Triplet(3 * T(i, j)+2, 3 * T(i, k)+1, m(3*j+2, 3*k+1)));
                triplets.push_back(Triplet(3 * T(i, j)+2, 3 * T(i, k)+2, m(3*j+2, 3*k+2)));
            }
        }
    }

    M.resize(3 * qdot.rows(), 3 * qdot.rows());
    M.setFromTriplets(triplets.begin(), triplets.end());
}
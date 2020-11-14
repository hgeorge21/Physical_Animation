#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    
    for (int i = 0; i < q.size(); i++) {
        triplets.push_back(T(i, i, mass));
    }

    M.resize(q.size(), q.size());
    M.setFromTriplets(triplets.begin(), triplets.end());
}

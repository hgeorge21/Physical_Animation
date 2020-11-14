#include <fixed_point_constraints.h>
#include <algorithm>
#include <numeric>

void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
	P.resize(q_size - 3 * indices.size(), q_size);

	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	std::vector<unsigned int> all_q(int(q_size / 3));
	std::iota(std::begin(all_q), std::end(all_q), 0);
	std::vector<unsigned int> non_fixed;

	std::set_difference(all_q.begin(), all_q.end(), indices.begin(), indices.end(), std::inserter(non_fixed, non_fixed.begin()));

	for (int i = 0; i < non_fixed.size(); i++) {
		triplets.push_back(T(3 * i, 3 * non_fixed[i], 1.0));
		triplets.push_back(T(3 * i + 1, 3 * non_fixed[i] + 1, 1.0));
		triplets.push_back(T(3 * i + 2, 3 * non_fixed[i] + 2, 1.0));
	}
	P.setFromTriplets(triplets.begin(), triplets.end());
}
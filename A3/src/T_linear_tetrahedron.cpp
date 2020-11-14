#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>

void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    Eigen::Matrix1212d M;
    mass_matrix_linear_tetrahedron(M, qdot, element, density, volume);

    Eigen::Vector12d qdot_;
    qdot_ << qdot.segment<3>(3 * element(0)), qdot.segment<3>(3 * element(1)), qdot.segment<3>(3 * element(2)), qdot.segment<3>(3 * element(3));

    T = 0.5 * qdot_.transpose() * M * qdot_;
}
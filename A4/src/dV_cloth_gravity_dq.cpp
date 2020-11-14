#include <dV_cloth_gravity_dq.h>

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    Eigen::VectorXd g_ = Eigen::VectorXd::Zero(fg.rows());
    for (int i = 0; i < M.rows() / 3; i++) {
        g_.segment<3>(3 * i) = g;
    }
    fg = -1. * M * g_;
}

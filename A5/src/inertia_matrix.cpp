#include <inertia_matrix.h>
#include <cassert>
#include <iostream>

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {
    // Compute center of mass
    double mass_x = 0;
    double mass_y = 0;
    center << 0., 0., 0.;
    for (int f = 0; f < F.rows(); f++) {
        Eigen::Vector3d X0 = V.row(F(f, 0));
        Eigen::Vector3d X1 = V.row(F(f, 1));
        Eigen::Vector3d X2 = V.row(F(f, 2));

        Eigen::Vector3d N = (X1 - X0).cross(X2 - X0).normalized();

        center(0) += (N(0) * (X0(0) * X0(0) + X0(0) * X1(0) + X0(0) * X2(0) + X1(0) * X1(0) + X1(0) * X2(0) + X2(0) * X2(0))) / 12;
        center(1) += (N(1) * (X0(1) * X0(1) + X0(1) * X1(1) + X0(1) * X2(1) + X1(1) * X1(1) + X1(1) * X2(1) + X2(1) * X2(1))) / 12;
        center(2) += (N(2) * (X0(2) * X0(2) + X0(2) * X1(2) + X0(2) * X2(2) + X1(2) * X1(2) + X1(2) * X2(2) + X2(2) * X2(2))) / 12;
        
        mass_x += (N(0) * (X0(0) + X1(0) + X2(0))) / 6;
        mass_y += (N(1) * (X0(1) + X1(1) + X2(1))) / 6;
    }
    center = density / 2 * center;
    mass_x = density * mass_x;
    mass_y = density * mass_y;
    
    mass = F.rows() * (center(0) - center(1)) / (mass_x - mass_y);
    center = 1. / mass * center;

    
}
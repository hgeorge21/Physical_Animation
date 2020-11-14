#include <rodrigues.h>
#include <cmath>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {
	R << 0, -exp(omega(2)), exp(omega(1)),
		exp(omega(2)), 0, -exp(omega(0)),
		-exp(omega(1)), exp(omega(0)), 0;
}
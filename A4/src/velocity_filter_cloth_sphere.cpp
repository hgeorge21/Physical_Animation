#include <velocity_filter_cloth_sphere.h>

void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {
    
    for (int idx = 0; idx < indices.size(); idx++) {
        int i = indices[idx];
        Eigen::Vector3d n = normals[idx];
        Eigen::Vector3d v = qdot.segment<3>(3 * i);

        double a = -1. * std::min(0., n.dot(v));
        qdot.segment<3>(3 * i) = v + a * n;
    }
}
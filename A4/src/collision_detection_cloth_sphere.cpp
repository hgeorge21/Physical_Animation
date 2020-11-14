#include <collision_detection_cloth_sphere.h>
#include <iostream>
void collision_detection_cloth_sphere(std::vector<unsigned int> &cloth_index, std::vector<Eigen::Vector3d> &normals, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Vector3d> center, double radius) {
    cloth_index.clear();
    normals.clear();

    int nv = q.rows() / 3;
    for (int i = 0; i < nv; i++) {
        Eigen::Vector3d x = q.segment<3>(3 * i);
        double d = (x - center).dot(x - center);

        if (d <= pow(radius, 2)) {
            cloth_index.push_back(i);
            
            Eigen::Vector3d n = (x - center).normalized();
            normals.push_back(n);
        }
    }
}
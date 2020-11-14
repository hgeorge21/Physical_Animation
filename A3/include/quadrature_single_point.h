#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template<typename Ret, typename Integrand_Func>
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q, 
                                               Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                               Integrand_Func integrand) {
    Eigen::Vector3d x0 = q.segment<3>(3 * element(0));
    Eigen::Vector3d x1 = q.segment<3>(3 * element(1));
    Eigen::Vector3d x2 = q.segment<3>(3 * element(2));
    Eigen::Vector3d x3 = q.segment<3>(3 * element(3));
    
    Eigen::Vector3d centroid = 0.25 * (x0 + x1 + x2 + x3);
    
    integrand(integrated, q, element, centroid);
    integrated = volume * integrated;
}


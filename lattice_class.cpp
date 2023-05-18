#include "lattice_class.h"
#include "helper.h"
#include <numeric>
#include <vector>
#include <iterator>


SquareLattice::FieldType SquareLattice::generate_random_field_value(){
    double inclination =   M_PI * ((double)rand())/((double)RAND_MAX);
    double azimuth =  2. * M_PI * ((double)rand())/((double)RAND_MAX);

    return {
        sin(inclination) * cos(azimuth),
        sin(inclination) * sin(azimuth),
        cos(inclination)
    };
}


// Sometimes this gives nan -> Why? Because of the way the phi are randomized?
SquareLattice::PhiType SquareLattice::compute_QL_for_single_triangle(SiteIndex s, int t) {
    FieldType phi_1 = phi[s];
    FieldType phi_2 = phi[triangles[s][t][0]];
    FieldType phi_3 = phi[triangles[s][t][1]];

    PhiType norm1 = sqrt(std::inner_product(phi_1.begin(), phi_1.end(), phi_1.begin(), 0.0));
    PhiType norm2 = sqrt(std::inner_product(phi_2.begin(), phi_2.end(), phi_2.begin(), 0.0));
    PhiType norm3 = sqrt(std::inner_product(phi_3.begin(), phi_3.end(), phi_3.begin(), 0.0));

    // Pre-compute dot products for re-use.
    PhiType dot_12 = std::inner_product(phi_1.begin(), phi_1.end(), phi_2.begin(), 0.0);
    PhiType dot_23 = std::inner_product(phi_2.begin(), phi_2.end(), phi_3.begin(), 0.0);
    PhiType dot_31 = std::inner_product(phi_3.begin(), phi_3.end(), phi_1.begin(), 0.0);

    PhiType rho = sqrt(2.0 * (1.+dot_12) * (1.+dot_23) * (1.+dot_31));
    PhiType ql_cos = (1. + dot_12 + dot_23 + dot_31) / rho;

    FieldType cross_product = outer_product(phi_2, phi_3);
    PhiType ql_sin = std::inner_product(phi_1.begin(), phi_1.end(), cross_product.begin(), 0.0)  / rho;

    // std:: cout << "----" << std::endl;
    // std::cout << "p1(" << norm1 << ") = [" << phi_1[0] << ", " << phi_1[1] << ", " << phi_1[2] << "]" << std::endl;
    // std::cout << "p2(" << norm2 << ")= [" << phi_2[0] << ", " << phi_2[1] << ", " << phi_2[2] << "]" << std::endl;
    // std::cout << "p3(" << norm3 << ") = [" << phi_3[0] << ", " << phi_3[1] << ", " << phi_3[2] << "]" << std::endl;
    // std::cout << dot_12 << " // " << dot_23 << " // " << dot_31 << std::endl;
    // std::cout << "Rho "<< rho << std::endl;
    // std::cout << "QL from acos: "<< acos(ql_cos) << " // " << ql_cos << std::endl;
    // std::cout << "QL from asin: "<< asin(ql_sin) << " // " << ql_sin << std::endl;
    // std::cout << "Diff: "<< asin(ql_sin) - acos(ql_cos) << " // " << ql_sin << std::endl;

    return acos(ql_cos) / (2*M_PI);
}


SquareLattice::PhiType SquareLattice:: compute_QL(SiteIndex s) {
  PhiType ql = 0.;
  for (int t=0; t<8; t++){
      ql += compute_QL_for_single_triangle(s, t);
  }
  return ql;
}

// Attention: There might be a factor 2 here because of double summation of the neighbors.
// This should be checked.
SquareLattice::PhiType SquareLattice:: compute_AL(SiteIndex s) {
  PhiType a_sum = 0.;
  for (int n=0; n<4; n++){
      a_sum += std::inner_product(phi[s].begin(), phi[s].end(), phi[neighbors[s][n]].begin(), 0.0);
      // cout << a_sum*beta << endl;
  }
  return -a_sum * beta;
}

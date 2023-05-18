#include "lattice_class.h"
#include "helper.h"
#include <iostream>
#include <numeric>
#include <vector>
#include <iterator>

void SquareLattice::init_lattice() {
    for (SiteIndex s=0; s < n_sites; s++) {
        phi.push_back(generate_random_unit_vector());
    }
}

// DEBUG STUFF.
void SquareLattice::dump() {
    for (SiteIndex s=0; s < n_sites; s++) {
//        double dot = std::inner_product(phi[s].begin(), phi[s].end(), phi[s].begin(), 0.0);
//        std::cout << dot << " [" << phi[s][0] << " // " << phi[s][1] << " // " << phi[s][2] << "]" << std::endl;
        std::cout << " [" << triangles[s][0][0] << " // " << triangles[s][0][1]  << "]" << std::endl;
    }
}

void SquareLattice::make_triangles(){
    for (int i=0; i<Lx; i++) {
        for (int j=0; j<Ly;j++){
            array<array<SiteIndex,2>,8> triangle;

            triangle[0] = {
                get_site(i, previous_y(j)),
                get_site(next_x(i), previous_y(j))
            };
            triangle[1] = {
                get_site(next_x(i), previous_y(j)),
                get_site(next_x(i), j)
            };
            triangle[2] = {
                get_site(next_x(i), j),
                get_site(next_x(i), next_y(j))
            };
            triangle[3] = {
                get_site(next_x(i), next_y(j)),
                get_site(i, next_y(j))
            };
            triangle[4] = {
                get_site(i, next_y(j)),
                get_site(previous_x(i), next_y(j))
            };
            triangle[5] = {
                get_site(previous_x(i), next_y(j)),
                get_site(previous_x(i), j)
            };
            triangle[6] = {
                get_site(previous_x(i), j),
                get_site(previous_x(i), previous_y(j))
            };
            triangle[7] = {
                get_site(previous_x(i), previous_y(j)),
                get_site(i, previous_y(j))
            };

            triangles.push_back(triangle);
        }
    }
}


SquareLattice::PhiType SquareLattice::compute_QL_for_single_triangle(SiteIndex s, int t) {
//    array<PhiType,3> phi_1 = phi[s];
//    array<PhiType,3> phi_2 = phi[triangles[s][t][0]];
//    array<PhiType,3> phi_3 = phi[triangles[s][t][1]];

    array<PhiType,3> phi_1 = generate_random_unit_vector();
    array<PhiType,3> phi_2 = generate_random_unit_vector();
    array<PhiType,3> phi_3 = generate_random_unit_vector();

//    array<PhiType,3> phi_1 = normalize({1, 2, 3});
//    array<PhiType,3> phi_2 = normalize({4, 5, 6});
//    array<PhiType,3> phi_3 = normalize({7, 8, 1});

    PhiType norm1 = sqrt(std::inner_product(phi_1.begin(), phi_1.end(), phi_1.begin(), 0.0));
    PhiType norm2 = sqrt(std::inner_product(phi_2.begin(), phi_2.end(), phi_2.begin(), 0.0));
    PhiType norm3 = sqrt(std::inner_product(phi_3.begin(), phi_3.end(), phi_3.begin(), 0.0));

    // Pre-compute dot products for re-use.
    PhiType dot_12 = std::inner_product(phi_1.begin(), phi_1.end(), phi_2.begin(), 0.0);
    PhiType dot_23 = std::inner_product(phi_2.begin(), phi_2.end(), phi_3.begin(), 0.0);
    PhiType dot_31 = std::inner_product(phi_3.begin(), phi_3.end(), phi_1.begin(), 0.0);

    PhiType rho = sqrt(2.0 * (1.+dot_12) * (1.+dot_23) * (1.+dot_31));
    PhiType ql_cos = (1. + dot_12 + dot_23 + dot_31) / rho;

    array<PhiType,3> cross_product = outer_product(phi_2, phi_3);
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

    return acos(ql_cos);
}


SquareLattice::PhiType SquareLattice:: compute_QL() {
    PhiType ql = 0.;

    for (SiteIndex s=0; s<n_sites; s++) {
        for (int t=0; t<8; t++){
            ql += compute_QL_for_single_triangle(s, t);
        }
    }
    return ql;
}

#include <vector>
#include <array>
#include <cmath>
#include <iostream>

using namespace std;

class SquareLattice {
    // Some type definitions for better readability.

public:
    using PhiType = double;
    using FieldType = array<PhiType,3>;
    using SiteIndex = int;

    SquareLattice(int Lx, int Ly, double beta, double theta): Lx(Lx), Ly(Ly), beta(beta), theta(theta) {
        n_sites = Lx * Ly;

        make_triangles();
        make_neighbours();

        init_random_fields();
    };

    // Used for diebugging.
    void dump() {
        for (SiteIndex s=0; s < n_sites; s++) {
    //        double dot = std::inner_product(phi[s].begin(), phi[s].end(), phi[s].begin(), 0.0);
    //        std::cout << dot << " [" << phi[s][0] << " // " << phi[s][1] << " // " << phi[s][2] << "]" << std::endl;
            std::cout << " [" << triangles[s][0][0] << " // " << triangles[s][0][1]  << "]" << std::endl;
        }
    }

    PhiType compute_action() { return compute_AL() + theta*compute_QL(); }

private:
    int Lx, Ly; // Extent of the lattice.
    double theta; // Topological term.
    double beta; // Inverse coupling.
    int n_sites; // Number of sites.
    vector<FieldType> phi; // Holds the values of phi.
    vector<array<array<SiteIndex,2>,8>> triangles; // Holds all 8 triangles for each site.
    vector<array<SiteIndex,4>> neighbors; // Holds all 4 neighbors for each site.

    SiteIndex next_x(int ix) { return next(ix, Lx); };
    SiteIndex next_y(int iy) { return next(iy, Ly); };
    SiteIndex next(int i, int L) {
        if (i == L) return 0;
        return i + 1;
    };

    SiteIndex previous_x(int ix) { return previous(ix, Lx); };
    SiteIndex previous_y(int iy) { return previous(iy, Ly); };
    SiteIndex previous(int i, int L) {
        if (i == 0) return L-1;
        return i - 1;
    };

    SiteIndex get_site(int i, int j) {
        return i + j*Lx;
    }

    void make_neighbours() {
        for (int i=0; i<Lx; i++) {
            for (int j=0; j<Ly; j++) {
                // std::cout << get_site(i, previous_y(j)) << " //";
                // std::cout << get_site(next_x(i), j) << " //";
                // std::cout << get_site(i, next_y(j)) << " //";
                // std::cout << get_site(previous_x(i), j) << " //" << std::endl;

                neighbors.push_back({
                    get_site(i, previous_y(j)),
                    get_site(next_x(i), j),
                    get_site(i, next_y(j)),
                    get_site(previous_x(i), j)
                });
            }
        }
    };

    void make_triangles() {
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
    };


    // ---------------

    FieldType generate_random_field_value();
    void init_random_fields() {
        for (SiteIndex s=0; s < n_sites; s++) {
            phi.push_back(generate_random_field_value());
        }
    };

    PhiType compute_QL(SiteIndex s);
    PhiType compute_QL_for_single_triangle(SiteIndex s, int t);
    PhiType compute_QL() {
      PhiType ql = 0.;
      for (SiteIndex s=0; s<n_sites; s++) {
        ql += compute_QL(s);
      }
      return ql;
    }

    PhiType compute_AL(SiteIndex s);
    PhiType compute_AL() {
      PhiType al = 0.;
      for (SiteIndex s=0; s<n_sites; s++) {
        al += compute_AL(s);
      }
      return al;
    };
};

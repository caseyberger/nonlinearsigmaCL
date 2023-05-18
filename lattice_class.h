#include <vector>
#include <array>
#include <cmath>

using namespace std;

class SquareLattice {
    // Some type definitions for better readability.

public:
    using PhiType = double;
    using SiteIndex = int;

    SquareLattice(int Lx, int Ly): Lx(Lx), Ly(Ly) {
        n_sites = Lx * Ly;
        init_lattice();
        make_triangles();

        PhiType ql = compute_QL();
//        PhiType ql = compute_QL_for_single_triangle(0, 0);
    };
    void dump();

private:
    int Lx, Ly; // Extent of the lattice.
    int n_sites; // Number of sites.
    vector<array<PhiType,3>> phi; // Holds the values of phi.
    vector<array<array<SiteIndex,2>,8>> triangles; // Holds all 8 triangles for each node.

    void init_lattice();
    void make_triangles();

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

    array<PhiType,3> generate_random_unit_vector(){
        double inclination =   M_PI * ((double)rand())/((double)RAND_MAX);
        double azimuth =  2. * M_PI * ((double)rand())/((double)RAND_MAX);

        return {
            sin(inclination) * cos(azimuth),
            sin(inclination) * sin(azimuth),
            cos(inclination)
        };
    }

    PhiType compute_QL_for_single_triangle(SiteIndex s, int t);
    PhiType compute_QL();
};

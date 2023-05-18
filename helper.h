#include <vector>
#include <array>
#include <numeric>

array<double,3> outer_product(array<double,3> v1, array<double,3> v2) {
    return {
        v1[1]*v2[2] - v1[2]*v2[1],
        -v1[0]*v2[2] + v1[2]*v2[0],
        v1[0]*v2[1] - v1[1]*v2[0],
    };
};

array<double,3> normalize(array<double,3> v) {
    double dot = sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
    return { v[0]/dot, v[1]/dot, v[2]/dot };
}
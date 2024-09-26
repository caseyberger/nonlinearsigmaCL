//Casey Berger
//Created: May 29, 2023
//Last edited: June 16, 2023
//all useful math operations stored here

#include <array>
#include "mathlib.h"

double dot(std::array<double, 3> a, std::array<double, 3> b){
    double s = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return s;
}

std::array<double, 3> cross(std::array<double, 3> a, std::array<double, 3> b){
    std::array<double, 3> c;
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
}
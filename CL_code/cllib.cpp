//Casey Berger
//Created: April 2, 2024
//Last edited: April 2, 2024
//Misc tools for performing complex Langevin calculations here

#include <array>
#include "mathlib.h"

 std::array<double, 2> complex_dot(std::array<double, 6> a, std::array<double, 6> b){
    //even indices are the real components, odd indices are the imaginary components
    std::array<double, 3> a_Re = {a[0], a[2], a[4]};
    std::array<double, 3> b_Re = {b[0], b[2], b[4]};
    std::array<double, 3> a_Im = {a[1], a[3], a[5]};
    std::array<double, 3> b_Im = {b[1], b[3], b[5]};
    std::array<double, 2> s = {0.,0.};
    //x components
    
    //y components
    
    //z components
    
    
    return s;
}

std::array<double, 6> complex_cross(std::array<double, 6> a, std::array<double, 6> b){
    std::array<double, 6> c = {0,0,0,0,0,0};
    return c;
}
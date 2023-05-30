//Casey Berger
//Created: May 29, 2023
//Edited: May 29, 2023
//all useful math operations stored here

#include "mathlib.h"

double dot(double a[3], double b[3]){
    double s = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return s;
}

double* cross(double a[3], double b[3]){
    static double c[3];
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
}
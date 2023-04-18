// Casey Berger
// Created: Mar 9, 2023
// Last edited: Apr 18, 2023

#include <string> //string
#include <iostream> //cout
#include <cmath> //M_PI, sin, cos, pow
#include <iomanip> //std::setw

#include "test_suite.h"
#include "lattice.h"
#include "action_suite.h"

/* 
Anything needed to test things while debugging -- 
print functions, testing values, etc

void check_phi_magnitude(double *** Lattice, int len); 
double phi_tot(double *** Lattice, int len); 
void print_lattice(double *** Lattice, int len);
void print_value(double *** Lattice, int i, int j, int len, double value);
void test_triangles(int i, int j, int len);
void test_QL(double QLcos, double QLsin);
*/

void check_phi_magnitude(double *** Lattice, int len)
{
    double phi_magnitude = 0.0;
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            phi_magnitude = pow(Lattice[i][j][0],2) + pow(Lattice[i][j][1],2) + pow(Lattice[i][j][2],2);
            if (phi_magnitude != 1.){
                std::cout << "Phi is not a unit vector. Magnitude = ";
                std::cout << std::setw(5) << phi_magnitude << "at location";
                std::cout << "i,j = " << i << "," << j << std::endl;
            }
        }
    }
}

double phi_tot(double *** Lattice, int len)
{
     //std::cout << "phi_tot" << std::endl;
    double phi = 0.0;
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            phi += pow(Lattice[i][j][0],2) + pow(Lattice[i][j][1],2) + pow(Lattice[i][j][2],2);
        }
    }
    return phi;
}

void print_lattice(double *** Lattice, int len)
{
    //std::cout << "print_lattice" << std::endl;
    //prints lattice to screen
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            std::cout << std::setw(5) << Lattice[i][j][0];
            std::cout << std::setw(5) << Lattice[i][j][1];
            std::cout << std::setw(5) << Lattice[i][j][2] << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


void print_value(double *** Lattice, int i, int j, int len, double value, std::string valname)
{
    //std::cout << "function: print_value in test_suite.cpp" << std::endl;
    //prints value calculated on lattice to screen
    std::cout << std::setw(2) << "i" << std::setw(2) << "j";
    std::cout << std::setw(16) << valname << std::endl;
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            std::cout << std::setw(2) << i << std::setw(2) << j;
            std::cout << std::setw(16) << value<< std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void test_triangles(int i, int j, int len)
{
    //testing triangles and plus one minus one
    int triangles[8][3][2];
    make_triangles(i,j,len,triangles);
    for (int n = 0; n<8;n++)
    {
        std::cout << std::setw(2) << i << std::setw(2) << j;
        std::cout << std::setw(2) << triangles[n][0][0] << std::setw(2) << triangles[n][0][1];
        std::cout << std::setw(2) << triangles[n][1][0] << std::setw(2) << triangles[n][1][1];
        std::cout << std::setw(2) << triangles[n][2][0] << std::setw(2) << triangles[n][2][1]<< std::endl;
    }
    std::cout << std::endl;
}

void test_QL(double *** Lattice, int i, int j, int len)
{
    int triangles[8][3][2];
    double QLsin, QLcos;
    make_triangles(i,j,len,triangles);
    for (int n = 0; n < 8; n++)
    {
        QLsin = QL_triangle(triangles[n], Lattice, true);
        QLsin = QL_triangle(triangles[n], Lattice, false);
        std::cout << std::setw(10) << "triangle " << std::setw(2) << n;
        std::cout << std::setw(10) << "QLcos = " << std::setw(10) << QLcos;
        std::cout << std::setw(10) << "QLsin = "<< std::setw(10) << QLsin << std::endl;
    }
}

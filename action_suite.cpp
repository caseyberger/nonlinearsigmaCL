// Casey Berger
// Created: Apr 18, 2023
// Last edited: Apr 18, 2023

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>

#include "action_suite.h"
#include "lattice.h"

/*
This file contains all the operations needed to calculate the 
lattice action, such as constructing the triangles, calculating
the topological charge Q_L and A_L
MAY NEED TO CALL LATTICE.H???
*/

void make_triangles(int i, int j, int len, int (&triangles)[8][3][2]){
    //returns the 8 triangles formed by the plaquettes surrounding the point you're on
    //you need to do this with nearest neighbors! Do you need to do a plus_one, minus_one function??
    //std::cout << "make_triangles" << std::endl;
    
    //triangle 1 
    triangles[0][0][0] = i;
    triangles[0][0][1] = j;
    triangles[0][1][0] = plus_one(i,len);
    triangles[0][1][1] = plus_one(j,len);
    triangles[0][2][0] = plus_one(i,len);
    triangles[0][2][1] = j;
    
    //triangle 2
    triangles[1][0][0] = i;
    triangles[1][0][1] = j;
    triangles[1][1][0] = i;
    triangles[1][1][1] = plus_one(j,len);
    triangles[1][2][0] = plus_one(i,len);
    triangles[1][2][1] = plus_one(j,len);
    
    //triangle 3
    triangles[2][0][0] = i;
    triangles[2][0][1] = j;
    triangles[2][1][0] = minus_one(i,len);
    triangles[2][1][1] = plus_one(j,len);
    triangles[2][2][0] = i;
    triangles[2][2][1] = plus_one(j,len);
    
    //triangle 4
    triangles[3][0][0] = i;
    triangles[3][0][1] = j;
    triangles[3][1][0] = minus_one(i,len);
    triangles[3][1][1] = j;
    triangles[3][2][0] = minus_one(i,len);
    triangles[3][2][1] = plus_one(j,len);
    
    //triangle 5
    triangles[4][0][0] = i;
    triangles[4][0][1] = j;
    triangles[4][1][0] = minus_one(i,len);
    triangles[4][1][1] = minus_one(j,len);
    triangles[4][2][0] = minus_one(i,len);
    triangles[4][2][1] = j;
    
    //triangle 6
    triangles[5][0][0] = i;
    triangles[5][0][1] = j;
    triangles[5][1][0] = i;
    triangles[5][1][1] = minus_one(j,len);
    triangles[5][2][0] = minus_one(i,len);
    triangles[5][2][1] = minus_one(j,len);
    
    //triangle 7
    triangles[6][0][0] = i;
    triangles[6][0][1] = j;
    triangles[6][1][0] = plus_one(i,len);
    triangles[6][1][1] = minus_one(j,len);
    triangles[6][2][0] = i;
    triangles[6][2][1] = minus_one(j,len);
    
    //triangle 8
    triangles[7][0][0] = i;
    triangles[7][0][1] = j;
    triangles[7][1][0] = plus_one(i,len);
    triangles[7][1][1] = j;
    triangles[7][2][0] = plus_one(i,len);
    triangles[7][2][1] = minus_one(j,len);
}

double QL_triangle(int current_triangle[3][2], double *** Lattice){
    //Calculates QL on a single triangle
    // note -- this is not returning the same value of QL from sin and cos... is there another way to get rid of the imaginary part?
    double phi2crossphi3[3];
    double rho, rho2, QLcos, QLsin;
    int i1,j1,i2,j2,i3,j3;
    i1 = current_triangle[0][0];
    j1 = current_triangle[0][1];
    i2 = current_triangle[1][0];
    j2 = current_triangle[1][1];
    i3 = current_triangle[2][0];
    j3 = current_triangle[2][1];
    rho2 = 2.*(1. + dot_product(Lattice[i1][j1], Lattice[i2][j2]))*(1. + dot_product(Lattice[i2][j2], Lattice[i3][j3]))*(1. + dot_product(Lattice[i3][j3], Lattice[i1][j1]));
    rho = sqrt(rho2);
    QLcos = (1. + dot_product(Lattice[i1][j1], Lattice[i2][j2]) + dot_product(Lattice[i2][j2], Lattice[i3][j3]) + dot_product(Lattice[i3][j3], Lattice[i1][j1]))/rho;
    cross_product(Lattice[i2][j2],Lattice[i3][j3],phi2crossphi3);
    QLsin = dot_product(Lattice[i1][j1],phi2crossphi3)/rho;
    std::cout << "QL_triangle"<< std::endl;
    //test_QL(acos(QLcos), asin(QLsin));
    return acos(QLcos);
}

/*
double A_lattice(double beta, double *** Lattice){
    //calculates the standard lattice action A_L
    double A_L = 0.0;
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            int i_nn,j_nn;
            //neighbor i+1,j
            i_nn = plus_one(i,len);
            j_nn = j;
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j_nn]);
            
            //neighbor i-1,j
            i_nn = minus_one(i,len);
            j_nn = j;
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j_nn]);
            
            //neighbor i,j+1
            i_nn = i;
            j_nn = plus_one(j,len);
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j_nn]);
            
            //neighbor i,j-1
            i_nn = i;
            j_nn = minus_one(j,len);
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j_nn]);
            }
        }
    return -1.*beta*A_L;
}
*/
/*
double Q_lattice(double *** Lattice){
    //calculates topological charge
    
    //How to do this...
    //You need to loop over 8 triangles
    //You will need a function to identify each triangle, maybe return all of them... an array?
    //Then you calculate rho squared (double) by summing over the triangles
    //Then you calculate the not renormalized Q from rho squared (double) and return it. You can renormalize it later with Z
    double Q_L = 0.0;
    int triangles[8][3][2];
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            make_triangles(i,j,len,triangles);
            for (int n = 0; n < 8; n++)
            {
                double QL_tri = QL_triangle(triangles[n], Lattice);
            }
        }
    }
    return Q_L;
}
*/
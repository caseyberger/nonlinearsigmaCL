// Casey Berger
// Created: Apr 18, 2023
// Last edited: Apr 25, 2023

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

void pick_phi(int i, int j, double (&phi)[3], double *** Lattice, bool old_lattice){
    if(old_lattice==true){
        //set phi
        phi[0] = Lattice[i][j][0];
        phi[1] = Lattice[i][j][1];
        phi[2] = Lattice[i][j][2];
    }
    else{
        //set phi
        phi[0] = Lattice[i][j][3];
        phi[1] = Lattice[i][j][4];
        phi[2] = Lattice[i][j][5];
    }
}

void make_triangles(int i, int j, int len, int (&triangles)[8][3][2]){
    //returns the 8 triangles formed by the plaquettes surrounding the point you're on
    //std::cout << "make_triangles" << std::endl;
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: make_triangles"<< std::endl;
#endif  
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

double QL_triangle(int current_triangle[3][2], double *** Lattice, bool arcsin, bool old_lattice){
    //Calculates QL on a single triangle
    // note -- this is not returning the same value of QL from sin and cos... is there another way to get rid of the imaginary part?
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: QL_triangle"<< std::endl;
#endif  
    double phi2crossphi3[3], phi1[3], phi2[3], phi3[3];
    double rho, rho2, QLcos, QLsin;
    int i1,j1,i2,j2,i3,j3;
    i1 = current_triangle[0][0];
    j1 = current_triangle[0][1];
    i2 = current_triangle[1][0];
    j2 = current_triangle[1][1];
    i3 = current_triangle[2][0];
    j3 = current_triangle[2][1];
    pick_phi(i1,j1, phi1, Lattice, old_lattice)//set phi1
    pick_phi(i2,j2, phi2, Lattice, old_lattice)//set phi2
    pick_phi(i3,j3, phi3, Lattice, old_lattice)//set phi3
        
    rho2 = 2.*(1. + dot_product(phi1, phi2))*(1. + dot_product(phi2, phi3))*(1. + dot_product(phi3, phi1));
    rho = sqrt(rho2);
    QLcos = (1. + dot_product(phi1, phi2) + dot_product(phi2, phi3) + dot_product(phi3, phi1))/rho;
    cross_product(phi2, phi3, phi2crossphi3);
    QLsin = dot_product(phi1,phi2crossphi3)/rho;
    if (arcsin){
        return asin(QLsin);
    }
    else{
        return acos(QLcos);
    }   
}

double Q_lattice(double *** Lattice, int len, bool old_lattice){
    //calculates topological charge
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: Q_lattice in action_suite"<< std::endl;
#endif  
    double Q_L = 0.0;
    int triangles[8][3][2];
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            make_triangles(i,j,len,triangles);
            for (int n = 0; n < 8; n++) //loop over 8 triangles
            {
                double QL_tri = QL_triangle(triangles[n], Lattice, old_lattice);
                Q_L += QL_tri;
            }
        }
    }
    return Q_L;//This Q_L is not renormalized. You can renormalize it later with Z
}

double A_lattice(double beta, double *** Lattice, int len, bool old_lattice){
    //calculates the standard lattice action A_L
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: A_lattice in action_suite"<< std::endl;
#endif  
    double A_L = 0.0;
    double phi[3], phi_nn[3];
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            int i_nn,j_nn;
            //neighbor i+1,j
            i_nn = plus_one(i,len);
            j_nn = j;
            pick_phi(i,j, phi, Lattice, old_lattice)//set phi
            pick_phi(i_nn,j_nn, phi_nn, Lattice, old_lattice)//set phi_nn
            A_L += dot_product(phi,phi_nn);
            
            /* only positive nn?
            //neighbor i-1,j
            i_nn = minus_one(i,len);
            j_nn = j;
            pick_phi(i,j, phi, Lattice, old_lattice)//set phi
            pick_phi(i_nn,j_nn, phi_nn, Lattice, old_lattice)//set phi_nn
            phi_nn[2] = Lattice[i_nn][jj_nn][2];
            A_L += dot_product(phi,phi_nn);
            */
            
            //neighbor i,j+1
            i_nn = i;
            j_nn = plus_one(j,len);
            pick_phi(i,j, phi, Lattice, old_lattice)//set phi
            pick_phi(i_nn,j_nn, phi_nn, Lattice, old_lattice)//set phi_nn
            A_L += dot_product(phi,phi_nn);
            /* only positive nn?
            //neighbor i,j-1
            i_nn = i;
            j_nn = minus_one(j,len);//set phi
            pick_phi(i,j, phi, Lattice, old_lattice)//set phi
            pick_phi(i_nn,j_nn, phi_nn, Lattice, old_lattice)//set phi_nn
            A_L += dot_product(phi,phi_nn);
            */
            }
        }
    return -1.*beta*A_L;
}

double S_lattice(double beta, double *** Lattice, int len, double itheta, bool old_lattice){
    //calculates the full lattice action S_L = A_L - i theta Q_L
    //not sure yet how to deal with the imaginary part, so right now I'm making one variable called itheta that will be real and analytically continued to imaginary values
#ifdef EXTREME_TESTING_MODE
    std:: cout << "Function: S_lattice in action_suite"<< std::endl;
    std::cout << "Note: not sure yet how to deal with the imaginary part, so right now I'm making one variable called itheta that will be real and analytically continued to imaginary values" << std::endl;
#endif  
    double S_L = 0.0;
    S_L = A_lattice(beta,Lattice,len, old_lattice);
    S_L += -1. * itheta * Q_lattice(Lattice, len, old_lattice);
    return S_L;
}
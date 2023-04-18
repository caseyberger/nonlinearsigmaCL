//#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include "action_suite.h"

/*
This file contains all the operations needed to calculate the 
lattice action, such as constructing the triangles, calculating
the topological charge Q_L and A_L
MAY NEED TO CALL LATTICE.H???
*/

void make_triangles(int i, int j, int len, int (&triangles)[8][3][2]){
    //returns the 8 triangles formed by the plaquettes surrounding the point you're on
    //you need to do this with nearest neighbors! Do you need to do a plus_one, minus_one function??
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
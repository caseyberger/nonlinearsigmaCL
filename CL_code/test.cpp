// Casey Berger
// Created: April 10, 2024
// Last edited: April 10, 2024
//
// takes input file. Run with ./test_nlsigma inputs
//

// include files
#include <time.h>//time (used to set random number seed, and to calculate dt)
#include <iostream> //cout
#include <cmath> //M_PI
#include <cstdlib> //rand
#include <fstream> //fout
#include <string> //string
#include <sstream> //stringstream for logfile
#include <omp.h> //openMP

//custom header files
#include "CL.h"
#include "cllib.h"
#include "mathlib.h"

using namespace std;
using complexlangevin::CL;

int main (int argc, char *argv[])
{
    cout << "Running testing suite for CL_simulation." << endl;
    cout << "Starting clock." << endl;
    time_t begin, end;
    double dt;
    srand(1723); //seed random number
    time(&begin);
    
     //Initalize the lattice - dynamically allocate the memory for the lattice
    cout << "Constructing lattice" << endl;    
    CL Lattice(argc, argv);//construct lattice
    int length = Lattice.GetLength();
    cout << "Lattice has length " << length << endl;
    cout << "Setting new lattice length to " << 2*length << endl;
    Lattice.SetLength(2*length);
    length = Lattice.GetLength();
    cout << "Lattice now has length " << length << endl;
   
    time(&end);
    dt = end - begin;
    cout << "Total time elapsed: " << dt/60. << " minutes." << endl;

    return 0;
}


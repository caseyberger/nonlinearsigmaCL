// Casey Berger
// Created: April 2, 2024
// Last edited: April 2, 2024
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
#include "cllib.h"
#include "mathlib.h"

using namespace std;

int main (int argc, char *argv[])
{
#ifdef TESTING_MODE
    cout << "Testing mode ON." << endl;
#endif
    cout << "Starting clock." << endl;
    time_t begin, end;
    double dt;
    srand(1723); //seed random number
    time(&begin);
    
    
    
   
    time(&end);
    dt = end - begin;
    cout << "Total time elapsed: " << dt/60. << " minutes." << endl;

    return 0;
}


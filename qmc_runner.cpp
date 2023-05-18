#include <iostream>
#include <cmath>
#include "lattice_class.h"


using namespace std;

int main (int argc, char *argv[])
{
    srand(1723); //seed random number
    double beta = 1.6;

    SquareLattice lattice = SquareLattice(5, 5);

//    lattice.dump();
}
#include <iostream>
#include <cmath>
#include "lattice_class.h"


using namespace std;

int main (int argc, char *argv[])
{
    srand(1723); //seed random number
    double beta = 1.6;
    double theta = 0.0;
    int L = 5;


    SquareLattice lattice = SquareLattice(L, L, beta, theta);

    double s = lattice.compute_action();
    cout << "Action: " << s << endl;
}

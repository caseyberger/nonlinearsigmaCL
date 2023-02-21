// Casey Berger
// Created: Feb 21, 2023
// Last edited: Feb 21, 2023
//

// include files
#include <time.h>//time (used to set random number seed)
#include <iostream> //cout
#include <cmath> //M_PI, sin, cos
#include <iomanip> //setw
#include <cstdlib> //rand
//#include <fstream>
//#include <vector>
//#include <string>
//#include <stdio.h>
//#include <sstream>

using namespace std;

//global variables
const int len = 4; //length of lattice
const int l_end = len-1; //last spot on lattice
const int num = len*len;  //total number of lattice sites

//function declaration
void make_lattice(double (&Lattice)[len][len][3]);
void print_lattice(double Lattice[len][len][3]);
//void write_to_file();

int main ()
{
    srand(time(NULL)); //seed random number
    double Lattice[len][len][3]; //stores lattice configuration
    make_lattice(Lattice);
    print_lattice(Lattice);
    //write_to_file(Lattice);
    return 0;
}

void make_lattice(double (&Lattice)[len][len][3])
{
    //cout << "make_lattice" << endl;
    //generates random configuration
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            //generate a random polar and azimuthal angle
            double inclination =   M_PI * ((double)rand())/((double)RAND_MAX); //polar angle = inclination
            double azimuth =  2. * M_PI * ((double)rand())/((double)RAND_MAX); //azimuthal angle = azimuth
            //create unit spin vector components from angles
            Lattice[i][j][0] = sin(inclination) * cos(azimuth);
            Lattice[i][j][1] = sin(inclination) * sin(azimuth);
            Lattice[i][j][2] = cos(inclination);
        }
    }
}

void print_lattice(double Lattice[len][len][3])
{
    //cout << "print_lattice" << endl;
    //prints lattice to screen
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            cout << setw(5) << Lattice[i][j][0];
            cout << setw(5) << Lattice[i][j][1];
            cout << setw(5) << Lattice[i][j][2] << endl;
        }
        cout << endl;
    }
    cout << endl;
}
/*
void write_to_file(vector<double> &exact_E, vector<double> &mc_E, vector<double> &temp_vec)
{
    //cout << "write_to_file" << endl;
    //output both solutions to a .txt file to open in gnuplot
    string fname = "mc_ising_data.txt";
    ofstream fout; //output stream
    fout.open(fname.c_str(),ios::out);
    
    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout << setw(20) << "# E/N exact" << setw(20) << "# E/N mc" << setw(20) << "T/J" << endl;
    for (unsigned int n = 0; n<exact_E.size(); n++)
    {
        fout.setf(ios::fixed);
        fout << setw(20) << exact_E[n] << setw(20) << 1.0*mc_E[n]/(1.0*J*num) << setw(20) << temp_vec[n]/(abs(J)) << endl;
    }
    fout.close();
    exact_E.clear();
    mc_E.clear();
    temp_vec.clear();
}*/
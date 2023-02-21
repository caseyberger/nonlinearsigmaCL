// Casey Berger
// Created: Feb 21, 2023
// Last edited: Feb 21, 2023
//

// include files
#include <time.h>//time (used to set random number seed)
#include <iostream> //cout
#include <cmath> //M_PI, sin, cos, pow
#include <iomanip> //setw
#include <cstdlib> //rand
#include <fstream> //fout
//#include <vector>
#include <string> //string
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
void print_value(double Lattice[len][len][3], double value[len][len]);
void calculate_phi_magnitude(double Lattice[len][len][3], double (&phi_magnitude)[len][len]);
//void write_to_file(double Lattice[len][len][3]);

int main ()
{
    srand(time(NULL)); //seed random number
    double Lattice[len][len][3]; //stores lattice configuration
    double phi_mag[len][len]; //stores size of unit vector at each lattice site
    
    
    make_lattice(Lattice); //fills lattice with phi values
    //print_lattice(Lattice);
    for (int n = 0; n<10; n++){
        calculate_phi_magnitude(Lattice, phi_mag);
        print_value(Lattice, phi_mag);
        make_lattice(Lattice);
    }
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


void print_value(double Lattice[len][len][3], double value[len][len])
{
    //cout << "print_value" << endl;
    //prints value calculated on lattice to screen
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            cout << setw(2) << i << setw(2) << j;
            cout << setw(10) << value[i][j]<< endl;
        }
        cout << endl;
    }
    cout << endl;
}


void calculate_phi_magnitude(double Lattice[len][len][3], double (&phi_magnitude)[len][len])
{
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            phi_magnitude[i][j] += pow(Lattice[i][j][0],2) + pow(Lattice[i][j][1],2) + pow(Lattice[i][j][2],2);
        }
    }
    return 0;
}

/*
void write_to_file(double Lattice[len][len][3])
{
    //cout << "write_to_file" << endl;
    //output both solutions to a .txt file to open in gnuplot
    string fname = "nonlinearsigma_data.txt";
    ofstream fout; //output stream
    fout.open(fname.c_str(),ios::out);
    
    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout << setw(10) << "step" << setw(10) << "|phi|" << endl;
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
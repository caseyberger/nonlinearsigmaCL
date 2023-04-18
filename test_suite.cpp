// Casey Berger
// Created: Mar 9, 2023
// Last edited: Apr 18, 2023

#include <string> //string
#include <iostream> //cout
#include <cmath> //M_PI, sin, cos, pow
#include <iomanip> //setw

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
                cout << "Phi is not a unit vector. Magnitude = ";
                cout << setw(5) << phi_magnitude << "at location";
                cout << "i,j = " << i << "," << j << endl;
            }
        }
    }
}

double phi_tot(double *** Lattice, int len)
{
     //cout << "phi_tot" << endl;
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


void print_value(double *** Lattice, int i, int j, int len, double value)
{
    //cout << "print_value" << endl;
    //prints value calculated on lattice to screen
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            cout << setw(2) << i << setw(2) << j;
            cout << setw(10) << value<< endl;
        }
        cout << endl;
    }
    cout << endl;
}

void test_triangles(int i, int j, int len)
{
    //testing triangles and plus one minus one
    int triangles[8][3][2];
    make_triangles(i,j,len,triangles);
    for (int n = 0; n<8;n++)
    {
        cout << setw(2) << i << setw(2) << j;
        cout << setw(2) << triangles[n][0][0] << setw(2) << triangles[n][0][1];
        cout << setw(2) << triangles[n][1][0] << setw(2) << triangles[n][1][1];
        cout << setw(2) << triangles[n][2][0] << setw(2) << triangles[n][2][1]<< endl;
    }
    cout << endl;
}

void test_QL(double QLcos, double QLsin)
{
    cout << setw(10) << "QLcos = " << setw(10) << QLcos << setw(10) << "QLsin = "<< setw(10) << QLsin << endl;   
}

void read_in_inputs(int argc, char *argv[],int &len, int &num, double &beta)
{
    //read in parameters
    cout << "Reading in parameters from input file" << endl;
    string str, filename;
    int n_params = 2;
    string inputs [2] = {"L","beta"};//read in keywords for parameters
    if (argc != 2){ //exits if input file is not given
        cerr << "Usage: ./nonlinearsigma input.txt"<< endl << "Exiting program" << endl;
        exit(10);
    }
    else{
        ifstream input_file(argv[1]);
        if (!input_file.is_open()){
            cerr << "input file cannot be opened";
            exit(10);
        }
        else{
            int count = 0;
            cout << "Starting param search in file."<<endl;
            while (count < n_params) {
                while (getline(input_file, str)) {
                    //search for params in input
                    size_t found = str.find(inputs[count]);
                    size_t start;
                    if (found != string::npos) {
                        start = str.find_last_of(' ');
                        inputs[count] = str.substr(start + 1);
                        count++;
                    }
                    else{
                        //if your inputs file doesn't have that parameter listed 
                        cerr << "parameter "<< inputs[count] << " not in input file.";
                        exit(10);
                    }
                }
            }
            len = stod(inputs[0]);
            num = len*len;
            beta = stod(inputs[1]);
            cout << "parameters acquired" <<endl;    
        }
    }
}
// Casey Berger
// Created: Feb 21, 2023
// Last edited: Mar 9, 2023
//
// takes input file. Run with ./nonlinearsigma inputs
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

//custom header files
#include "test_suite.h"
#include "lattice.h"
#include "action_suite.h"

using namespace std;

//function declaration
double Z_renorm(double beta, int len);
void create_logfile();
void write_to_file(int n, double phi, double Q_L, double A_L, double S_L);
void read_in_inputs(int argc, char *argv[],int &len, int &num, double &beta);

int main (int argc, char *argv[])
{
    int len, num;
    double beta = 1.6;
    
    read_in_inputs(argc, argv,len, num, beta);
    cout << "len = " << len << endl;
    cout << "beta = " << beta << endl;
    
    //Initalize the lattice - dynamically allocate the memory for the lattice
    cout << "Allocating memory for lattice" << endl;
    double *** Lattice = new double**[num];
    for(int i = 0; i < len; i++){
        Lattice[i] = new double*[len];
    }
    //allocation - 3 phi components
    //note: can this be put into a separate file as well? Think about how that would work
    for(int i = 0; i < len; i++){
        for (int j = 0; j<len; j++){
            Lattice[i][j] = new double[3];
        }
    }
    lattice_init(Lattice, len);//initialize phi everywhere
        
    double phi_mag[len][len]; //stores size of unit vector at each lattice site
    
    create_logfile(); //generates logfile with header
    //print_lattice(Lattice, len);
    
    double phi = 0.0;
    double A_L = 0.0;
    double Q_L = 0.0;
    double S_L = 0.0;
    double itheta = M_PI;
    
    for (int n = 0; n<10; n++){
        phi = phi_tot(Lattice, len);
        A_L = A_lattice(beta, Lattice, len);
        Q_L = Q_lattice(Lattice, len);
        S_L = S_lattice(beta, Lattice, len, itheta);
        write_to_file(n, phi, Q_L, A_L, S_L);
        lattice_init(Lattice, len);
    }
    
    return 0;
}


double Z_renorm(double beta, int len){
    //pulls renormalization factor from table given in paper
    if (beta == 1.5 and len==120){
        return 0.285;
    }
    else if (beta == 1.6 and len==180){
        return 0.325;
    }
    else if (beta == 1.7 and len==340){
        return 0.380;
    }
    else if (beta == 1.75 and len==470){
        return 0.412;
    }
    else{
        return 0.;
    }
}

void create_logfile()
{
    //cout << "create_logfile" << endl;
    //create header of logfile before
    string fname = "nonlinearsigma_data.csv";
    ofstream fout; //output stream
    fout.open(fname.c_str(),ios::out);
    
    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout  << "step" << "," << "|phi|" << "," << "Q_L"<< "," << "A_L"<< "," << "S_L"<< endl;
    fout.close();
}

void write_to_file(int n, double phi, double Q_L, double A_L, double S_L)
{
    //cout << "write_to_file" << endl;
    //output calculations .csv file
    string fname = "nonlinearsigma_data.csv";
    ofstream fout; //output stream
    fout.open(fname.c_str(),std::ios_base::app);

    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout << n << "," << phi<< "," << Q_L<< "," << A_L <<"," << S_L << endl;
    fout.close();
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
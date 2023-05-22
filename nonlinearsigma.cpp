// Casey Berger
// Created: Feb 21, 2023
// Last edited: May 17, 2023
//
// takes input file. Run with ./nonlinearsigma inputs
//

// include files
#include <time.h>//time (used to set random number seed, and to calculate dt)
#include <iostream> //cout
#include <cmath> //M_PI, sin, cos, pow
#include <iomanip> //setw
#include <cstdlib> //rand
#include <fstream> //fout
#include <string> //string
//#include <vector>
//#include <stdio.h>
//#include <sstream>

//custom header files
#include "test_suite.h"
#include "lattice.h"
#include "action_suite.h"
#include "monte_carlo.h"

using namespace std;

//function declaration
double Z_renorm(double beta, int len);
void create_logfile();
void write_to_file(double dt, int n, double phi, double Q_L, double A_L, double S_L);
void read_in_inputs(int argc, char *argv[],int &len, int &num, int &ntherm, int &nMC, double &beta);

int main (int argc, char *argv[])
{
#ifdef TESTING_MODE
    cout << "Testing mode ON." << endl;
    cout << "Starting clock." << endl;
#endif
    srand(1723); //seed random number
    time_t begin, end, begin_therm, end_therm, begin_mc, dt_end, dt_start, end_mc;
    double dt;
    
    time(&begin);
    
#ifdef EXTREME_TESTING_MODE
    cout << "Testing mode is EXTREME." << endl;
#endif

    int len, num, ntherm, nMC;
    double beta = 1.6;
    
    read_in_inputs(argc, argv,len, num, ntherm, nMC, beta);
    cout << "len = " << len << endl;
    cout << "beta = " << beta << endl;
    
    //Initalize the lattice - dynamically allocate the memory for the lattice
#ifdef TESTING_MODE
    cout << "Allocating memory for lattice" << endl;
#endif
    double *** Lattice = new double**[num];
    for(int i = 0; i < len; i++){
        Lattice[i] = new double*[len];
    }
    //allocation - 3 phi components x 2 (old and new)
    for(int i = 0; i < len; i++){
        for (int j = 0; j<len; j++){
            Lattice[i][j] = new double[6];
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
    bool old_lattice = true;
    
    //thermalization loop
#ifdef TESTING_MODE
    cout << "Starting thermalization loop of length " << ntherm << endl;
#endif
    
    time(&begin_therm);
    for (int n = 0; n<ntherm; n++){
        //some sort of updating function in here
        Metropolis_loop(beta, itheta, Lattice, len);
    }
    time(&end_therm);
    
    dt = end_therm - begin_therm;
    //MC loop
#ifdef TESTING_MODE
    cout << "Thermalization loop duration: " << dt << " seconds."<< endl;
    cout << "Starting Monte Carlo loop of length " << nMC << endl;
#endif
    
    time(&begin_mc);

    for (int n = 0; n<nMC; n++){
        //some sort of updating function in here
        time(&dt_start);
        phi = phi_tot(Lattice, len, old_lattice);
        A_L = A_lattice(beta, Lattice, len, old_lattice);
        Q_L = Q_lattice(Lattice, len, old_lattice);
        S_L = S_lattice(beta, Lattice, len, itheta, old_lattice);
        time(&dt_end);
        dt = dt_end-dt_start;
        write_to_file(dt, n, phi, Q_L, A_L, S_L);
        lattice_init(Lattice, len);
    }
    time(&end_mc);
    
    dt = end_mc - begin_mc;
#ifdef TESTING_MODE
     cout << "MC loop duration: " << dt << " seconds." << endl;
#endif
    time(&end);
    dt = end - begin;
#ifdef TESTING_MODE
     cout << "Total time elapsed: " << dt << " seconds." << endl;
#endif
    return 0;
}


double Z_renorm(double beta, int len){
    //pulls renormalization factor from table given in paper
#ifdef TESTING_MODE
    cout << "Function: Z_renorm" << endl;
#endif
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
#ifdef TESTING_MODE
    cout << "Function: create_logfile" << endl;
#endif
    //create header of logfile
    string fname = "nonlinearsigma_data.csv";
#ifdef TESTING_MODE
    cout << "logfile name: " << fname << endl;
#endif
    ofstream fout; //output stream
    fout.open(fname.c_str(),ios::out);
    
    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout << "dt"<< "," << "step" << "," << "|phi|" << "," << "Q_L"<< "," << "A_L"<< "," << "S_L"<< endl;
    fout.close();
}

void write_to_file(double dt, int n, double phi, double Q_L, double A_L, double S_L)
{
#ifdef TESTING_MODE
    cout << "Function: write_to_file" << endl;
#endif
    //output calculations .csv file
    string fname = "nonlinearsigma_data.csv";
#ifdef TESTING_MODE
    cout << "Filename: " << fname << endl;
#endif
    ofstream fout; //output stream
    fout.open(fname.c_str(),std::ios_base::app);

    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout << dt <<","<< n << "," << phi<< "," << Q_L<< "," << A_L <<"," << S_L << endl;
    fout.close();
}

void read_in_inputs(int argc, char *argv[],int &len, int &num, int &ntherm, int &nMC, double &beta)
{
    //read in parameters
#ifdef TESTING_MODE
    cout << "Function: read_in_inputs" << endl;
#endif
    string str, filename;
    int n_params = 2;
    string inputs [4] = {"L","beta", "ntherm","nMC"};//read in keywords for parameters
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
#ifdef TESTING_MODE
    cout << "Starting param search in file: ";
    for (int n=0; n<n_params; n++){
        cout << inputs[n] << ',';
    }
    cout << endl;
#endif  
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
            ntherm = stod(inputs[2]);
            nMC = stod(inputs[3]);
#ifdef TESTING_MODE
    cout << "parameters acquired: ";
    for (int n=0; n<n_params; n++){
        cout << inputs[n] << ',';
    }
    cout << endl;
#endif  
        }
    }
}
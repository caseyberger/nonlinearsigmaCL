// Casey Berger
// Created: Feb 21, 2023
// Last edited: May 30, 2023
//
// takes input file. Run with ./nonlinearsigma inputs
//

// include files
#include <time.h>//time (used to set random number seed, and to calculate dt)
#include <iostream> //cout
#include <cmath> //M_PI
#include <cstdlib> //rand
#include <fstream> //fout
#include <string> //string
//#include <vector>
//#include <stdio.h>
//#include <sstream>

//custom header files
#include "lattice.h"

using namespace std;
using nonlinearsigma::Lattice;

//function declaration
double Z_renorm(double beta, int len);
void create_logfile();
void write_to_file(double dt, int n, double phi, double Q_L, double A_L, double S_L, double acc);
void read_in_inputs(int argc, char *argv[],int &len, int &num, int &ntherm, int &nMC, double &beta,double &itheta);
void test_phi_distribution(Lattice L);

int main (int argc, char *argv[])
{
#ifdef TESTING_MODE
    cout << "Testing mode ON." << endl;
    cout << "Starting clock." << endl;
#endif
    time_t begin, end, begin_therm, end_therm, begin_mc, dt_end, dt_start, end_mc;
    double dt;
    srand(1723); //seed random number
    time(&begin);
    
#ifdef EXTREME_TESTING_MODE
    cout << "Testing mode is EXTREME." << endl;
#endif

    int len, num, ntherm, nMC;
    double beta = 1.6;
    double itheta = M_PI;
    
    read_in_inputs(argc, argv,len, num, ntherm, nMC, beta, itheta);
    cout << "len = " << len << endl;
    cout << "beta = " << beta << endl;
    cout << "itheta = " << itheta << endl;
    
    //Initalize the lattice - dynamically allocate the memory for the lattice
#ifdef TESTING_MODE
    cout << "Constructing lattice" << endl;
#endif
    
    Lattice L(len, beta, itheta);//construct lattice
    
#ifdef TESTING_MODE
    cout << "Length = " << L.getLength() << endl;
    cout << "Resetting length" << endl;
    L.setLength(10);
    cout << "New length = " << L.getLength() << endl;
    
    cout << "Beta = " << L.getBeta() << endl;
    cout << "Resetting beta" << endl;
    L.setBeta(1.0);
    cout << "New beta = " << L.getBeta() << endl;
    
    cout << "itheta = " << L.getiTheta() << endl;
    cout << "Resetting itheta" << endl;
    L.setiTheta(0.5*M_PI);
    cout << "New itheta = " << L.getiTheta() << endl;
    //note -- you must reinitialize the lattice after changing length
#endif
    
#ifdef TESTING_MODE
    cout << "Initializing lattice" << endl;
#endif
    L.initialize(); //initialize 3-component phi everywhere
    
    create_logfile(); //generates logfile with header
    

#ifdef TESTING_MODE
    cout << "Printing lattice" <<endl;
    L.printLattice();

    cout << "Printing triangles" <<endl;
    for(int i = 0; i < len; i++){
        for (int j = 0; j<len; j++){
            L.printTriangles(i,j);
        }
    }
#endif   
    
    double phi = 0.0;
    double A_L = 0.0;
    double Q_L = 0.0;
    double S_L = 0.0;
    double acc = 0.0;
    
    //thermalization loop
#ifdef TESTING_MODE
    cout << "Starting thermalization loop of length " << ntherm << endl;
#endif
    
    time(&begin_therm);
    L.zeroCount();
    L.thermalize(ntherm);
    time(&end_therm);
    
    dt = end_therm - begin_therm;
    //MC loop
#ifdef TESTING_MODE
    cout << "Thermalization loop duration: " << dt << " seconds."<< endl;
    cout << "Starting Monte Carlo loop of length " << nMC << endl;
#endif
    
    time(&begin_mc);

    for (int n = 0; n<nMC; n++){
        time(&dt_start);
        L.metropolisStep();
#ifdef TESTING_MODE
        L.printLattice();
        test_phi_distribution(L);
#endif
        phi = L.getPhiTot();
        A_L = L.calcAL();
        Q_L = L.calcQL();
        S_L = L.calcSL();
        acc = L.acceptanceRate();
        time(&dt_end);
        dt = dt_end-dt_start;
        write_to_file(dt, n, phi, Q_L, A_L, S_L, acc);
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
    fout << "dt,step,|phi|,Q_L,A_L,S_L,acc" << endl;
    fout.close();
}

void write_to_file(double dt, int n, double phi, double Q_L, double A_L, double S_L, double acc)
{
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
    fout << dt <<","<< n << "," << phi<< "," << Q_L<< "," << A_L << "," << S_L << "," << acc << endl;
    fout.close();
}

void read_in_inputs(int argc, char *argv[],int &len, int &num, int &ntherm, int &nMC, double &beta, double &itheta)
{
    //read in parameters -- note itheta is read in as a multiple/fraction of pi
#ifdef TESTING_MODE
    cout << "Function: read_in_inputs" << endl;
#endif
    string str, filename;
    const int n_params = 5;
    string inputs [n_params] = {"L","beta", "itheta", "ntherm","nMC"};//read in keywords for parameters
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
            itheta = stod(inputs[2])*M_PI;
            ntherm = stod(inputs[3]);
            nMC = stod(inputs[4]);
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

void test_phi_distribution(Lattice L){
    //output phi distributions as .csv file
    cout << "Saving phi distribution to file" <<endl;
    //create header of logfile
    string fname = "phi_test.csv";
    ofstream fout; //output stream
    fout.open(fname.c_str(),ios::out);
    
    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout << "i,j,phi_x,phi_y,phi_z,r1,r2" << endl;
    int len = 100; //in order to get a large sample -- maybe make larger?
    L.setLength(len);
    L.initialize();
    for (int i = 0; i < len; i++){
        for (int j = 0; j < len; j++){
            double *phi = L.getPhi(i,j);
            double *r = L.getRandNums();
            fout << i <<","<< j << "," << phi[0]<< "," << phi[1]<< "," << phi[2] << "," << r[0] << "," << r[1] << endl;
        }
    }
    fout.close();
}
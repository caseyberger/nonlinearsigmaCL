// Casey Berger
// Created: Feb 21, 2023
// Last edited: June 12, 2023
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
#include <sstream> //stringstream for logfile

//custom header files
#include "lattice.h"

using namespace std;
using nonlinearsigma::Lattice;

//function declaration
double Z_renorm(double beta, int len);
void create_logfile(Lattice L);
void write_to_file(Lattice L, int n, double phi, double Q_L, double A_L, double S_L, double Xi_L, double F_L[2], double acc);
void read_in_inputs(int argc, char *argv[],int &len, int &ntherm, int &nMC, double &beta,double &itheta);
void test_phi_distribution(Lattice L);
void save_correlation_function(Lattice L);
void testing_suite(int len, double beta, double itheta);

int main (int argc, char *argv[])
{
#ifdef TESTING_MODE
    cout << "Testing mode ON." << endl;
#endif
    cout << "Starting clock." << endl;
    time_t begin, end, begin_therm, end_therm, begin_mc, time_now, end_mc;
    double dt;
    srand(1723); //seed random number
    time(&begin);
    
    int len = 180;
    int ntherm = 1000;
    int nMC = 100;
    double beta = 1.6;
    double itheta = M_PI;
    
    read_in_inputs(argc, argv,len, ntherm, nMC, beta, itheta);
    cout << "len = " << len << endl;
    cout << "beta = " << beta << endl;
    cout << "itheta = " << itheta << endl;
    cout << "ntherm = " << ntherm << endl;
    cout << "nMC = " << nMC << endl;
     
#ifdef EXTREME_TESTING_MODE
    testing_suite(len, beta, itheta);
    cout << "Testing concluded" << endl;
    exit(0);
#endif
    
    //Initalize the lattice - dynamically allocate the memory for the lattice
    cout << "Constructing lattice" << endl;    
    Lattice L(len, beta, itheta);//construct lattice
    L.setnTherm(ntherm);
    L.setnMC(nMC);

    cout << "Initializing lattice" << endl;
    L.initialize(); //initialize 3-component phi everywhere
    create_logfile(L); //generates logfile with header 
    
    double phi = 0.0;
    double A_L = 0.0;
    double Q_L = 0.0;
    double S_L = 0.0;
    double Xi_L = 0.0;
    double* F_L;
    double acc = 0.0;
    
    //thermalization loop
    cout << "Starting thermalization loop of length " << L.getnTherm() << endl;
    
    time(&begin_therm);
    L.thermalize();
    time(&end_therm);
    
    dt = end_therm - begin_therm;
    
    //MC loop
    cout << "Thermalization loop duration: " << dt/60. << " minutes."<< endl;
    cout << "Starting Monte Carlo loop of length " << L.getnMC() << endl;
    
    time(&begin_mc);

    for (int n = 0; n<L.getnMC(); n++){
        L.metropolisStep();
        phi  = L.getPhiTot();
        A_L  = L.calcAL();
        Q_L  = L.calcQL();
        S_L  = L.calcSL();
        Xi_L = L.calcXi();
        F_L  = L.calcF();
        acc  = L.acceptanceRate();
        write_to_file(L, n, phi, Q_L, A_L, S_L, Xi_L, F_L, acc);
#ifdef TESTING_MODE
        time(&time_now);
        if (n%100 == 0){
            cout << "On MC step " << n << ", total time elapsed = " << time_now - begin << " seconds." << endl;
        }
#endif
    }
    save_correlation_function(L);
    time(&end_mc);
    
    dt = end_mc - begin_mc;
    cout << "MC loop duration: " << dt/60. << " minutes." << endl;

    time(&end);
    dt = end - begin;
    cout << "Total time elapsed: " << dt/60. << " minutes." << endl;

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


void create_logfile(Lattice L)
{
    string fname = L.getFilename();
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
    fout << "step,|phi|,Q_L,A_L,S_L,Xi_L,F_LRe, F_LIm, acc" << endl;
    fout.close();
}

void write_to_file(Lattice L, int n, double phi, double Q_L, double A_L, double S_L, double Xi_L, double F_L[2], double acc)
{
    //output calculations .csv file
    string fname = L.getFilename();
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
    fout << n << "," << phi<< "," << Q_L<< "," << A_L << ",";
    fout << S_L << ","<< Xi_L << "," << F_L[0] << "," << F_L[1] << "," << acc << endl;
    fout.close();
}

void read_in_inputs(int argc, char *argv[],int &len, int &ntherm, int &nMC, double &beta, double &itheta)
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
    cout << "Saving phi magnitude and distribution to file" <<endl;
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
    fout << "i,j,phi_x,phi_y,phi_z,|phi|,r1,r2" << endl;
    int len = 100; //in order to get a large sample -- maybe make larger?
    L.setLength(len);
    L.initialize();
    for (int i = 0; i < len; i++){
        for (int j = 0; j < len; j++){
            field phi = L.getPhi(i,j);
            double phimag = L.getPhiMag(i,j);
            double *r = L.getRandNums();
            fout << i <<","<< j << ",";
            fout << phi[0]<< "," << phi[1]<< "," << phi[2] << "," << phimag << ",";
            fout << r[0] << "," << r[1] << endl;
        }
    }
    fout.close();
}

void save_correlation_function(Lattice L){
    //output phi distributions as .csv file
    cout << "Saving average correlation function to file" <<endl;
    //create header of logfile
    string fname = "Gij_avg_"+L.getFilename();
    ofstream fout; //output stream
    fout.open(fname.c_str(),ios::out);
    
    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout << "i,j,G_avg" << endl;
    int len = L.getLength(); 
    for (int i = 0; i < len; i++){
        for (int j = 0; j < len; j++){
            double G = L.getAvgG(i,j);
            fout << i <<","<< j << ",";
            fout << G << endl;
        }
    }
    fout.close();
}
    
void testing_suite(int len, double beta, double itheta){
    cout << "Extreme testing mode enabled. Running through tests." << endl;
    
    //construct lattice
    cout << "Constructing lattice" << endl;
    Lattice L(len, beta, itheta);
    
    //test setting and getting parameters
    cout << "Length = " << L.getLength() << endl;
    cout << "Resetting length" << endl;
    L.setLength(len + 1);
    cout << "New length = " << L.getLength() << endl;
    len = L.getLength();
    
    cout << "Beta = " << L.getBeta() << endl;
    cout << "Resetting beta" << endl;
    L.setBeta(1.0);
    cout << "New beta = " << L.getBeta() << endl;
    beta = L.getBeta();
    
    cout << "itheta = " << L.getiTheta() << endl;
    cout << "Resetting itheta" << endl;
    L.setiTheta(0.5*M_PI);
    cout << "New itheta = " << L.getiTheta() << endl;
    itheta = L.getiTheta();
    
    cout << "nTherm = " << L.getnTherm() << endl;
    cout << "Resetting nTherm" << endl;
    L.setnTherm(10);
    cout << "New nTherm = " << L.getnTherm() << endl;
    int ntherm = L.getnTherm();
    
    cout << "nMC = " << L.getnMC() << endl;
    cout << "Resetting nMC" << endl;
    L.setnMC(100);
    cout << "New nMC = " << L.getnMC() << endl;
    int nMC = L.getnMC();
    
    //initialization tests
    cout << "Initializing lattice" << endl;
    L.initialize();
    cout << endl;
    
    //print lattice
    cout << "Printing lattice" << endl;
    L.printLattice();
    
    //testing distrubution of phi
    test_phi_distribution(L);
    
    cout << "Total phi on the lattice is " << L.getPhiTot() << ", and total lattice sites is " << len*len << endl << endl;
    
    cout << "Resetting length" << endl;
    L.setLength(6);
    cout << "New length = " << L.getLength() << endl;
    len = L.getLength();
    
    //testing triangles
    cout << "Testing triangle generation" <<endl;
    int i = 0;
    int j = 0;
    L.printTriangles(i,j);
    i = len/2;
    L.printTriangles(i,j);
    j = len/2;
    L.printTriangles(i,j);
    cout << endl;
    
    //testing neighbor getting functions
    cout << "Testing ability to get neighboring phis" << endl;
    int testlen = 3;
    L.setLength(testlen);
    L.printLattice();
    for (int i = 0; i < testlen; i++){
        for (int j = 0; j < testlen; j++){
            L.printNeighbors(i,j);
        }
    }
    cout << endl;
    //testing lattice quantities
    cout << "Testing lattice calculations. " << endl;
    double QL, AL, SL;
    for(int i = 0; i < len; i++){
        for (int j = 0; j<len; j++){
            L.checkQL(i, j);
        }
    }
    L.setLength(2);
    L.setBeta(1.);
    L.setiTheta(M_PI);
    L.initialize();
    L.printLattice();
    AL = L.calcAL();
    QL = L.calcQL();
    SL = L.calcSL();
    cout << "for L = 2, beta = 1, itheta = pi, and phi pointing in random direction:" << endl;
    cout << "AL = " << AL << ", QL = " << QL << ", SL = AL - itheta QL = " << SL << endl << endl;
    
    L.setiTheta(0.);
    L.initialize();
    L.printLattice();
    AL = L.calcAL();
    QL = L.calcQL();
    SL = L.calcSL();
    cout << "for L = 2, beta = 1, itheta = 0, and phi pointing in random direction:" << endl;
    cout << "AL = " << AL << ", QL = " << QL << ", SL = AL - itheta QL = " << SL << endl << endl;
    
    L.setiTheta(M_PI);
    L.fixRNG(0.,0.);
    L.initialize();
    L.printLattice();
    AL = L.calcAL();
    QL = L.calcQL();
    SL = L.calcSL();
    cout << "for L = 2, beta = 1, itheta = pi, and all phi pointing in +z direction:" << endl;
    cout << "AL = " << AL << ", QL = " << QL << ", SL = AL - itheta QL = " << SL << endl << endl;
   
    L.fixRNG(0.5,0.);
    L.initialize();
    L.printLattice();
    AL = L.calcAL();
    QL = L.calcQL();
    SL = L.calcSL();
    cout << "for L = 2, beta = 1, itheta = pi, and all phi pointing in +x direction:" << endl;
    cout << "AL = " << AL << ", QL = " << QL << ", SL = AL - itheta QL = " << SL << endl << endl;
    
    L.fixRNG(0.5,0.5);
    L.initialize();
    L.printLattice();
    AL = L.calcAL();
    QL = L.calcQL();
    SL = L.calcSL();
    cout << "for L = 2, beta = 1, itheta = pi, and all phi pointing in -x direction:" << endl;
    cout << "AL = " << AL << ", QL = " << QL << ", SL = AL - itheta QL = " << SL << endl << endl;
    
    L.fixRNG(0.5,0.25);
    L.initialize();
    L.printLattice();
    AL = L.calcAL();
    QL = L.calcQL();
    SL = L.calcSL();
    cout << "for L = 2, beta = 1, itheta = pi, and all phi pointing in y direction:" << endl;
    cout << "AL = " << AL << ", QL = " << QL << ", SL = AL - itheta QL = " << SL << endl << endl;
    
    L.fixRNG(0.5,0.75);
    L.initialize();
    L.printLattice();
    AL = L.calcAL();
    QL = L.calcQL();
    SL = L.calcSL();
    cout << "for L = 2, beta = 1, itheta = pi, and all phi pointing in -y direction:" << endl;
    cout << "AL = " << AL << ", QL = " << QL << ", SL = AL - itheta QL = " << SL << endl << endl;
    
    L.freeRNG();
    L.initialize();
    L.printLattice();
    AL = L.calcAL();
    QL = L.calcQL();
    SL = L.calcSL();
    cout << "for L = 2, beta = 1, itheta = pi, and phi pointing in random direction:" << endl;
    cout << "AL = " << AL << ", QL = " << QL << ", SL = AL - itheta QL = " << SL << endl << endl;
    
    
    //MC testing
    L.setLength(4);
    L.setBeta(1.6);
    L.setiTheta(0.);
    L.setnTherm(10);
    L.initialize();
    L.printLattice();
    cout << "Metropolis testing, step 1" << endl;
    L.metropolisStep();
    L.printLattice();
    cout << "Metropolis testing, step 2" << endl;
    L.metropolisStep();
    L.printLattice();
    cout << "Metropolis testing, step 3" << endl;
    L.metropolisStep();
    L.printLattice();
    cout << "Metropolis testing, step 4" << endl;
    L.metropolisStep();
    L.printLattice();
    cout << "Acceptance rate = " << L.acceptanceRate() << endl;
    L.zeroCount();
    L.thermalize();
    cout << "Acceptance rate = " << L.acceptanceRate() << endl;
}
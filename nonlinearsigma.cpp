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
//#include "test_suite.h"
#include "lattice.h"

using namespace std;

//global variables
//const int len = 4; //length of lattice
//const int num = len*len;  //total number of lattice sites

//function declaration
void check_phi_magnitude(double *** Lattice, int len);
double phi_tot(double *** Lattice, int len); //only useful for testing
double dot_product(double vec1[3], double vec2[3]);
void cross_product(double vec1[3], double vec2[3], double (&cross_prod)[3]);
int plus_one(int i, int len);
int minus_one(int i, int len);
//double A_lattice(double beta, double *** Lattice);
void make_triangles(int i, int j, int len, int (&triangles)[8][3][2]);
//double QL_triangle(int current_triangle[3][2], double *** Lattice);
//double Q_lattice(double *** Lattice);
double Z_renorm(double beta, int len);
void create_logfile();
void write_to_file(int n, double phi, double A_L);
void print_lattice(double *** Lattice, int len);
void print_value(double *** Lattice, int i, int j, int len, double value);
void test_triangles(int i, int j, int len);
void test_QL(double QLcos, double QLsin);
void read_in_inputs(int argc, char *argv[],int &len, int &num, double &beta);

int main (int argc, char *argv[])
{

    
    int len, num;
    double beta = 1.6;
    
    read_in_inputs(argc, argv,len, num, beta);
    
    //Initalize the lattice - dynamically allocate the memory for the lattice
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
    
    for (int n = 0; n<10; n++){
        phi = phi_tot(Lattice, len);
        //A_L = A_lattice(beta, Lattice);
        //Q_L = Q_lattice(Lattice);
        write_to_file(n, phi, A_L);
        lattice_init(Lattice, len);
    }
    
    return 0;
}

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

double dot_product(double vec1[3], double vec2[3]){
    //calculates the dot product of two vectors
    double dot_prod = 0.0;
    dot_prod = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
    return dot_prod;
}

void cross_product(double vec1[3], double vec2[3],double (&cross_prod)[3]){
    //calculates the cross product of two vectors
    cross_prod[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    cross_prod[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
    cross_prod[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

int plus_one(int i, int len){
    //returns site plus one, using periodic boundary conditions
    if (i==len-1){
        return 0;
    }
    else{
        return i+1;
    }
}

int minus_one(int i, int len){
    //returns site minus one, using periodic boundary conditions
    if (i==0){
        return len-1;
    }
    else{
        return i-1;
    }
}
/*
double A_lattice(double beta, double *** Lattice){
    //calculates the standard lattice action A_L
    double A_L = 0.0;
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            int i_nn,j_nn;
            //neighbor i+1,j
            i_nn = plus_one(i,len);
            j_nn = j;
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j_nn]);
            
            //neighbor i-1,j
            i_nn = minus_one(i,len);
            j_nn = j;
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j_nn]);
            
            //neighbor i,j+1
            i_nn = i;
            j_nn = plus_one(j,len);
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j_nn]);
            
            //neighbor i,j-1
            i_nn = i;
            j_nn = minus_one(j,len);
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j_nn]);
            }
        }
    return -1.*beta*A_L;
}
*/
void make_triangles(int i, int j, int len, int (&triangles)[8][3][2]){
    //returns the 8 triangles formed by the plaquettes surrounding the point you're on
    //you need to do this with nearest neighbors! Do you need to do a plus_one, minus_one function??
    //triangle 1 
    triangles[0][0][0] = i;
    triangles[0][0][1] = j;
    triangles[0][1][0] = plus_one(i,len);
    triangles[0][1][1] = plus_one(j,len);
    triangles[0][2][0] = plus_one(i,len);
    triangles[0][2][1] = j;
    
    //triangle 2
    triangles[1][0][0] = i;
    triangles[1][0][1] = j;
    triangles[1][1][0] = i;
    triangles[1][1][1] = plus_one(j,len);
    triangles[1][2][0] = plus_one(i,len);
    triangles[1][2][1] = plus_one(j,len);
    
    //triangle 3
    triangles[2][0][0] = i;
    triangles[2][0][1] = j;
    triangles[2][1][0] = minus_one(i,len);
    triangles[2][1][1] = plus_one(j,len);
    triangles[2][2][0] = i;
    triangles[2][2][1] = plus_one(j,len);
    
    //triangle 4
    triangles[3][0][0] = i;
    triangles[3][0][1] = j;
    triangles[3][1][0] = minus_one(i,len);
    triangles[3][1][1] = j;
    triangles[3][2][0] = minus_one(i,len);
    triangles[3][2][1] = plus_one(j,len);
    
    //triangle 5
    triangles[4][0][0] = i;
    triangles[4][0][1] = j;
    triangles[4][1][0] = minus_one(i,len);
    triangles[4][1][1] = minus_one(j,len);
    triangles[4][2][0] = minus_one(i,len);
    triangles[4][2][1] = j;
    
    //triangle 6
    triangles[5][0][0] = i;
    triangles[5][0][1] = j;
    triangles[5][1][0] = i;
    triangles[5][1][1] = minus_one(j,len);
    triangles[5][2][0] = minus_one(i,len);
    triangles[5][2][1] = minus_one(j,len);
    
    //triangle 7
    triangles[6][0][0] = i;
    triangles[6][0][1] = j;
    triangles[6][1][0] = plus_one(i,len);
    triangles[6][1][1] = minus_one(j,len);
    triangles[6][2][0] = i;
    triangles[6][2][1] = minus_one(j,len);
    
    //triangle 8
    triangles[7][0][0] = i;
    triangles[7][0][1] = j;
    triangles[7][1][0] = plus_one(i,len);
    triangles[7][1][1] = j;
    triangles[7][2][0] = plus_one(i,len);
    triangles[7][2][1] = minus_one(j,len);
}
/*
double QL_triangle(int current_triangle[3][2], double *** Lattice){
    //Calculates QL on a single triangle
    // note -- this is not returning the same value of QL from sin and cos... is there another way to get rid of the imaginary part?
    double phi2crossphi3[3];
    double rho, rho2, QLcos, QLsin;
    int i1,j1,i2,j2,i3,j3;
    i1 = current_triangle[0][0];
    j1 = current_triangle[0][1];
    i2 = current_triangle[1][0];
    j2 = current_triangle[1][1];
    i3 = current_triangle[2][0];
    j3 = current_triangle[2][1];
    rho2 = 2.*(1. + dot_product(Lattice[i1][j1], Lattice[i2][j2]))*(1. + dot_product(Lattice[i2][j2], Lattice[i3][j3]))*(1. + dot_product(Lattice[i3][j3], Lattice[i1][j1]));
    rho = sqrt(rho2);
    QLcos = (1. + dot_product(Lattice[i1][j1], Lattice[i2][j2]) + dot_product(Lattice[i2][j2], Lattice[i3][j3]) + dot_product(Lattice[i3][j3], Lattice[i1][j1]))/rho;
    cross_product(Lattice[i2][j2],Lattice[i3][j3],phi2crossphi3);
    QLsin = dot_product(Lattice[i1][j1],phi2crossphi3)/rho;
    cout << "QL_triangle"<< endl;
    test_QL(acos(QLcos), asin(QLsin));
    return 0;
}
*/
/*
double Q_lattice(double *** Lattice){
    //calculates topological charge
    
    //How to do this...
    //You need to loop over 8 triangles
    //You will need a function to identify each triangle, maybe return all of them... an array?
    //Then you calculate rho squared (double) by summing over the triangles
    //Then you calculate the not renormalized Q from rho squared (double) and return it. You can renormalize it later with Z
    double Q_L = 0.0;
    int triangles[8][3][2];
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            make_triangles(i,j,len,triangles);
            for (int n = 0; n < 8; n++)
            {
                double QL_tri = QL_triangle(triangles[n], Lattice);
            }
        }
    }
    return Q_L;
}*/
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
    fout << setw(10) << "step" << setw(10) << "|phi|" << setw(10) << "A_L"<< endl;
    fout.close();
}

void write_to_file(int n, double phi, double A_L)
{
    //cout << "write_to_file" << endl;
    //output both solutions to a .txt file to open in gnuplot
    string fname = "nonlinearsigma_data.txt";
    ofstream fout; //output stream
    fout.open(fname.c_str(),std::ios_base::app);

    // check if files are open
    if (!fout.is_open())
    {
        cerr << "Unable to open file " << fname <<"." << endl;
        exit(10);
    }
    fout.setf(ios::fixed);
    fout << setw(10) << n << setw(10) << phi<< setw(10) << A_L << endl;
    fout.close();
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
    string str, filename;
    int n_params = 2;
    string inputs [2] = {"L","beta"};//read in keywords for parameters
    if (argc != 2){ //exits if input file is not given
        cerr << "Usage: ./v1.4 input_CLB.txt"<<endl << "Exiting program" << endl;
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
                }
            }
            len = stod(inputs[0]);
            num = len*len;
            beta = stod(inputs[1]);
            //cout << "parameters acquired" <<endl;    
        }
    }
}
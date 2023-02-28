// Casey Berger
// Created: Feb 21, 2023
// Last edited: Feb 28, 2023
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
void calculate_phi_magnitude(double Lattice[len][len][3], double (&phi_magnitude)[len][len]);//consider changing to a "check" function that returns an error?
double phi_tot(double Lattice[len][len][3]); //only useful for testing
double dot_product(double vec1[3], double vec2[3]);
void find_nn(int (&nn_arr)[2]);
double A_lattice(double gL, double Lattice[len][len][3]);
void create_logfile();
void write_to_file(int n, double phi);

int main ()
{
    srand(time(NULL)); //seed random number
    double Lattice[len][len][3]; //stores lattice configuration
    double phi_mag[len][len]; //stores size of unit vector at each lattice site
    
    make_lattice(Lattice); //fills lattice with phi values
    create_logfile(); //generates logfile with header
    //print_lattice(Lattice);
    
    double phi = 0.0;
    double gL = 1./1.6;
    double A_L = 0.0;
    
    for (int n = 0; n<10; n++){
        phi = phi_tot(Lattice);
        A_L = A_lattice(gL, Lattice)
        write_to_file(n, phi);
        make_lattice(Lattice);
    }
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
            phi_magnitude[i][j] = pow(Lattice[i][j][0],2) + pow(Lattice[i][j][1],2) + pow(Lattice[i][j][2],2);
        }
    }
}

double phi_tot(double Lattice[len][len][3])
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

void find_nn(int (&nn_arr)[2]){
    //returns the phi nearest neighbor specified
    int i_nn = nn_arr[0];
    int j_nn = nn_arr[1];
    if (j_nn==len){
        j_nn = 0;
    }
    else if (j_nn ==-1){
        j_nn = len-1;
    }
    if (i_nn == len){
        i_nn = 0;
    }
    else if(i_nn == -1){
        i_nn = len-1;
    }
    nn_arr[0] = i_nn;
    nn_arr[1] = j_nn;
}

double A_lattice(double gL, double Lattice[len][len][3]){
    //calculates the standard lattice action A_L
    double A_L = 0.0;
    for (int i = 0; i<len; i++)
    {
        for (int j = 0; j<len; j++)
        {
            int i_nn,j_nn;
            int nn_arr[2];
            //neighbor i+1,j
            nn_arr[0] = i+1;
            nn_arr[1] = j;            
            find_nn(nn_arr);
            i_nn = nn_arr[0];
            j_nn = nn_arr[1];
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j-nn]);
            
            //neighbor i-1,j
            nn_arr[0] = i-1;
            nn_arr[1] = j;            
            find_nn(nn_arr);
            i_nn = nn_arr[0];
            j_nn = nn_arr[1];
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j-nn]);
        
            //neighbor i,j+1
            nn_arr[0] = i;
            nn_arr[1] = j+1;            
            find_nn(nn_arr);
            i_nn = nn_arr[0];
            j_nn = nn_arr[1];
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j-nn]);
            
            //neighbor i,j-1
            nn_arr[0] = i;
            nn_arr[1] = j-1;            
            find_nn(nn_arr);
            i_nn = nn_arr[0];
            j_nn = nn_arr[1];
            A_L += dot_product(Lattice[i][j],Lattice[i_nn][j-nn]);
            }
        }
    return -1.*A_L/gL;
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
    fout << setw(10) << "step" << setw(10) << "|phi|" << endl;
    fout.close();
}

void write_to_file(int n, double phi)
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
    fout << setw(10) << n << setw(10) << phi << endl;
    fout.close();
}
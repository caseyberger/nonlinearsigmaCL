// Casey Berger
// Created: Mar 9, 2023
// Last edited: Mar 9, 2023

#include <string> //string
#include <iostream> //cout
#include <cmath> //M_PI, sin, cos, pow
#include <iomanip> //setw

//global variables
const int len = 4; //length of lattice
const int l_end = len-1; //last spot on lattice
const int num = len*len;  //total number of lattice sites

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

void test_triangles(int i, int j)
{
    //testing triangles and plus one minus one
    int triangles[8][3][2];
    make_triangles(i,j,triangles);
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
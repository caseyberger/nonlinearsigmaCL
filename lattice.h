// Casey Berger
// Created: May 24, 2023
// Last edited: May 30, 2023

#pragma once

namespace nonlinearsigma{
    class Lattice {
        public:
        //constructor
        Lattice(int length, double beta, double itheta);
        //other public functions
        //set or update lattice parameters
        void setLength(int length);
        void setBeta(double beta);
        void setiTheta(double itheta);
        
        //retrieve lattice parameters
        int getLength();
        double getBeta();
        double getiTheta();
        double* getPhi(int i, int j);
        double getPhiMag(int i, int j);
        double getPhiTot();
        
        //initialize the lattice
        void initialize(); //tested 5/30/2023
        void testLattice();
        
        //print things to the screen
        void printLattice();
        void printTriangles(int i, int j);
        
        //calculate lattice quantities
        double calcQL();
        double calcAL();
        double calcSL();
        
        //monte carlo tools
        void metropolisStep();
        void thermalize(int ntherm);
        void zeroCount();
        double acceptanceRate();
    
        //private members -- only accessible within the class functions
        private:
        //variables
        int length_;
        double ***grid_;
        int *****triangles_;
        double beta_;
        double itheta_;
        int acceptCount_;
        int rejectCount_;
        double accRate_;
        
        //functions
        double* makePhi_();
        int plusOne_(int i);
        int minusOne_(int i);
        void makeTriangles_();
        int*** trianglesCCW_(int i, int j);
        double locQL_(int i, int j, int n, bool use_arccos);
        int* getNeighbors_(int i, int j);
        double** getNeighborPhis_(int i, int j);
        void printPhi_(int i, int j);
    };  
}



// Casey Berger
// Created: May 24, 2023
// Last edited: May 29, 2023
// This needs to: define a phi at each point in the lattice, pick out phi values given a location on the lattice, calculate dot and cross products.
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
        
        //initialize the lattice
        void initialize();
        
        //calculating lattice quantities
        double calcQL();
        double calcAL();
        double calcSL();
        
    
        private:
        int length_;
        double ***grid;
        int *****triangles_;
        double beta_;
        double itheta_;
        double* makePhi_();
        int plusOne_(int i);
        int minusOne_(int i);
        void makeTriangles_();
        int*** trianglesCCW_(int i, int j);
        double locQL_(int i, int j, int n, bool use_arccos);
        int* getNeighbors_(int i, int j);
        int** getNeighborPhis_(int i, int j);
    };  
}



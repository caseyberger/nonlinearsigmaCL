// Casey Berger
// Created: May 24, 2023
// Last edited: June 7, 2023

#pragma once

namespace nonlinearsigma{
    class Lattice {
        public:
        //constructor
        Lattice(int length, double beta, double itheta);
        //other public functions
        //set or update lattice parameters
        void setLength(int length);//tested 6/1/2023
        void setBeta(double beta);//tested 6/1/2023
        void setiTheta(double itheta);//tested 6/1/2023
        void setPhi(int i, int j, double phi[3]);//tested 6/5/2023
        void setAvgG(int i, int j, double Gij);
        void setnTherm(int ntherm);
        void setnMC(int nMC);
        void fixRNG(double r1, double r2);//tested 6/1/2023
        void freeRNG();//tested 6/1/2023
        
        //retrieve lattice parameters
        int getLength();//tested 6/1/2023
        double getBeta();//tested 6/1/2023
        double getiTheta();//tested 6/1/2023
        int getnTherm();
        int getnMC();
        string getFilename();
        double* getPhi(int i, int j);//tested 5/30/2023
        double* getRandNums();//tested 6/1/2023
        double getPhiMag(int i, int j); //tested 6/1/2023
        double getPhiTot(); //tested 6/1/2023
        double getAvgG(int i, int j);
        
        //initialize the lattice
        void initialize(); //tested 5/30/2023
        
        //print things to the screen
        void printLattice();//tested 5/30/2023
        void printTriangles(int i, int j); //tested 6/1/2023
        void printNeighbors(int i, int j);//tested 6/1/2023
        
        //calculate and test lattice quantities
        double calcQL();//tested 6/1/2023 -- note it's not producing integers!!
        void checkQL(int i, int j);//tested 6/1/2023
        double calcAL();//tested 6/1/2023
        double calcSL();//tested 6/1/2023
        double twoPointG(int i, int j);
        double calcXi();
        double* calcF();
        
        //monte carlo tools
        void metropolisStep();//tested 6/5/2023
        void thermalize();//tested 6/5/2023
        void zeroCount();//tested 6/5/2023
        double acceptanceRate();//tested 6/5/2023
    
        //private members -- only accessible within the class functions
        private:
        //variables
        int length_;
        double ***grid_;
        int *****triangles_;
        double beta_;
        double itheta_;
        double r1_;
        double r2_;
        bool fixedr_;
        int nTherm_;
        int nMC_;
        int acceptCount_;
        int rejectCount_;
        double accRate_;
        double **Gij_;
        string filename_;
        
        //functions
        double* makePhi_(); //tested 6/1/2023
        int plusOne_(int i); //tested 5/30/2023
        int minusOne_(int i);//tested 5/30/2023
        void makeTriangles_(); //tested 6/1/2023
        int*** trianglesCCW_(int i, int j); //tested 6/1/2023
        double locQL_(int i, int j, int n, bool use_arccos);//tested 6/1/2023
        int* getNeighbors_(int i, int j);//tested 6/1/2023
        double** getNeighborPhis_(int i, int j);//tested 6/1/2023
        void printPhi_(int i, int j); //tested 5/30/2023
        void generateFilename_();
    };  
}



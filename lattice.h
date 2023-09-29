// Casey Berger
// Created: May 24, 2023
// Last edited: July 14, 2023 - triangles
#include <array>
#include <vector>
#include <omp.h>
#pragma once


namespace nonlinearsigma{
    class Lattice {
        public:
        
        using vertex = std::array<int,2>;
        using triangle = std::array<vertex,3>;
        using site_triangles = std::array<triangle,2>;
        using field = std::array<double, 3>;
        
        //constructor
        Lattice(int length, double beta, double itheta);
        //other public functions
        //set or update lattice parameters
        void setLength(int length);//tested 6/1/2023
        void setBeta(double beta);//tested 6/1/2023
        void setiTheta(double itheta);//tested 6/1/2023
        void setPhi(int i, int j, field phi);//tested 6/5/2023
        void setnTherm(int ntherm);
        void setnMC(int nMC);
        void setFreq(int freq);
        void setTrig(bool use_arcsin);
        void fixRNG(double r1, double r2);//tested 6/1/2023
        void freeRNG();//tested 6/1/2023
        
        //retrieve lattice parameters
        int getLength();//tested 6/1/2023
        double getBeta();//tested 6/1/2023
        double getiTheta();//tested 6/1/2023
        int getnTherm();
        int getnMC();
        std::string getFilename();
        field getPhi(int i, int j);//tested 5/30/2023
        double* getRandNums();//tested 6/1/2023
        double getPhiMag(int i, int j); //tested 6/1/2023
        double getPhiTot(); //tested 6/1/2023
        double getAvgG(int i, int j);
        
        //initialize the lattice
        void initialize(); //tested 5/30/2023
        
        //print things to the screen
        void printLattice();//tested 5/30/2023
        void printTriangles(int i, int j);//updated 7/14/2023
        void printNeighbors(int i, int j);//tested 6/1/2023
        
        //write things to file
        void saveConfig(int step);
        
        //calculate and test lattice quantities
        double calcQL();//updated 7/14/2023
        void checkQL();//updated 7/14/2023
        double calcAL();//tested 6/1/2023
        double calcSL();//tested 6/1/2023
        void calcGij();
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
        std::vector < std::vector < field > > grid_;
        std::vector < std::vector < site_triangles > > triangles_;
        double beta_;
        double itheta_;
        double r1_;
        double r2_;
        bool fixedr_;
        bool use_arcsin_;
        int nTherm_;
        int nMC_;
        int freq_;
        int acceptCount_;
        int rejectCount_;
        double accRate_;
        std::vector < std::vector < double > > Gij_;
        std::string filename_;
        
        //functions
        field makePhi_(); //tested 6/1/2023
        int plusOne_(int i); //tested 5/30/2023
        int minusOne_(int i);//tested 5/30/2023
        void makeTriangles_();//updated 7/14/2023
        double locQL_(int i, int j, int n, bool use_arccos);//updated 7/14/2023
        std::array < vertex, 4 > getNeighbors_(int i, int j);//tested 6/1/2023
        std::array < field, 4 > getNeighborPhis_(int i, int j);//tested 6/1/2023
        void printPhi_(int i, int j); //tested 5/30/2023
        void generateFilename_();
    };  
}
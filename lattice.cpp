// Casey Berger
// Created: Mar 28, 2023
// Last edited: May 30, 2023

#include <iostream> //cout, endl
#include <cmath> //sqrt, sin, cos, acos, asin, exp
#include "mathlib.h" //dot, cross

#include "lattice.h"


namespace nonlinearsigma{
    //public functions
    //constructor
    Lattice::Lattice(int length, double beta, double itheta){
        Lattice::setLength(length); //set length
        Lattice::setBeta(beta); //set beta
        Lattice::setiTheta(itheta); //set itheta
    }
    //other public functions
    void Lattice::setLength(int length){
        length_ = length;
    }
    
    void Lattice::setBeta(double beta){
        beta_ = beta;
    }
    
    void Lattice::setiTheta(double itheta){
        itheta_ = itheta;
    }
    
    int Lattice::getLength(){
        return length_;
    }
    
    double Lattice::getBeta(){
        return beta_;
    }
    
    double Lattice::getiTheta(){
        return itheta_;
    }
    
    double* Lattice::getPhi(int i, int j){
        return grid_[i][j];
    }
    
    double Lattice::getPhiMag(int i, int j){
        double *phi = Lattice::getPhi(i,j);
        double phi_mag = dot(phi,phi);
        return phi_mag;
    }
    
    double Lattice::getPhiTot(){
        double phi_tot = 0.;
        for(int i = 0; i < length_; i++){
            for(int j = 0; j < length_; j++){
                phi_tot += Lattice::getPhiMag(i, j);
            }
        }
        return phi_tot;
    }
    
    void Lattice::initialize(){
        double *** grid = new double**[length_];
        for(int i = 0; i < length_; i++){
            grid[i] = new double*[length_];
        }
        //allocation - 3 phi components
        for(int i = 0; i < length_; i++){
            for (int j = 0; j<length_; j++){
                grid[i][j] = new double[3];
                grid[i][j] = Lattice::makePhi_();
            }
        }
        grid_ = grid;
        
        Lattice::makeTriangles_();
    }
    
    void Lattice::printLattice(){
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                double *phi = Lattice::getPhi(i,j);
                std::cout << "At point (" << i << "," << j <<"), phi = (";
                std::cout << phi[0] << "," << phi[1] << ","<< phi[2] << ")" << std::endl;
            }
        }
    }
    
    void Lattice::printTriangles(int i, int j){
        std::cout << "At point (" << i << "," << j << ")," << std::endl;
        for (int n = 0; n < 8; n++){
            std::cout << "Triangle "<< n + 1 << " = (";
            std::cout << triangles_[i][j][n][0][0] << ","<< triangles_[i][j][n][0][1] << "), (";
            std::cout << triangles_[i][j][n][1][0] << ","<< triangles_[i][j][n][1][1] << "), (";
            std::cout << triangles_[i][j][n][2][0] << ","<< triangles_[i][j][n][2][1]<< ")" << std::endl;
        }
    }
    
    double Lattice::calcQL(){
        //calculates topological charge  
        double Q_L = 0.0;
        bool use_arccos = true;//uses arccos to find QL for each triangle
        for (int i = 0; i<length_; i++)
        {
            for (int j = 0; j<length_; j++)
            {
                for (int n = 0; n < 8; n++) //loop over 8 triangles
                {
                    Q_L += Lattice::locQL_(i, j, n, use_arccos);
                }
            }
        }
        return Q_L;//This Q_L is not renormalized. You can renormalize it later with Z
    }
    
    double Lattice::calcAL(){
        //calculates the standard lattice action A_L
        double A_L = 0.0;
        for (int i = 0; i<length_; i++)
        {
            for (int j = 0; j<length_; j++)
            {
                double *phi = Lattice::getPhi(i,j);
                double **phiNN = Lattice::getNeighborPhis_(i,j); //0 and 1 are + direction
                //nearest neighbors in positive direction:
                A_L += dot(phi, phiNN[0]) + dot(phi, phiNN[1]);
                //nearest neighbors in negative direction:
                A_L += dot(phi, phiNN[2]) + dot(phi, phiNN[2]);
            }
        }
        return -1.*beta_*A_L;
    }
    
    double Lattice::calcSL(){
        //calculates the full lattice action S_L = A_L - i theta Q_L
        //not sure yet how to deal with the imaginary part, 
        //so right now I'm making one variable called itheta that will be real 
        //and analytically continued to imaginary values 
        double S_L = Lattice::calcAL() - 1.*itheta_*Lattice::calcQL();
        return S_L;
    }
    
    void Lattice::metropolisStep(){
        double Si, Sf, dS, r;
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                Si = Lattice::calcSL();
                double *phi_old = Lattice::getPhi(i, j);
                //update lattice
                double *phi_new = Lattice::makePhi_();
                grid_[i][j] = phi_new;
                Sf = Lattice::calcSL();
                dS = Sf - Si;
#ifdef TEST_CONSTANT_RN
                r = 0.5;
#else
                r = ((double)std::rand())/((double)RAND_MAX);
#endif
                if(dS < 0 || std::exp(dS) > r){
                    acceptCount_++;//increment accept counter
#ifdef TESTING_MODE
                    std:: cout << "Accept" << std::endl;
#endif
                }
                else{
                    grid_[i][j] = phi_old;//change the value back to the old phi
                    rejectCount_++;//increment reject counter
#ifdef TESTING_MODE
                    std:: cout << "Reject" << std::endl;
#endif
                }
            }//loop over j
        }//loop over i
        double acc_rate = (double)acceptCount_/((double)acceptCount_ + (double)rejectCount_);
        accRate_ = acc_rate;
    }
    
    void Lattice::thermalize(int ntherm){
        for (int n = 0; n < ntherm; n++){
            Lattice::metropolisStep();
        }
        //should you zero the count at the end of this?
    }
    
    void Lattice::zeroCount(){
        acceptCount_ = 0;
        rejectCount_ = 0;
    }
    
    double Lattice::acceptanceRate(){
        return accRate_;
    }
    
    //private functions
    double* Lattice::makePhi_(){
        static double phi[3];
        double r1, r2;
#ifdef TEST_CONSTANT_RN
        r1 = 0.5;
        r2 = 0.5;
#else
        r1 = ((double)std::rand())/((double)RAND_MAX);
        r2 = ((double)std::rand())/((double)RAND_MAX);
#endif
        
#ifdef TESTING_MODE
        std::cout << "r1 = " << r1 << ", r2 = " << r2 << std::endl;
#endif
        //generate a random polar and azimuthal angle
        double inclination =   M_PI * r1; //polar angle = inclination
        double azimuth =  2. * M_PI * r2; //azimuthal angle = azimuth
        //create unit spin vector components from angles
        phi[0] = std::sin(inclination) * std::cos(azimuth);
        phi[1] = std::sin(inclination) * std::sin(azimuth);
        phi[2] = std::cos(inclination);
        
        return phi;
    }
    
    int Lattice::plusOne_(int i){
        //int len = Lattice::getLength();
        if(i+1 == length_){ return 0;}
        else{return i+1;}
    }
    
    int Lattice::minusOne_(int i){
        if(i == 0){ return length_-1;}
        else{return i-1;}
    }
    
    void Lattice::makeTriangles_(){
        int ***** triangles = new int****[length_];
        for(int i = 0; i < length_; i++){
            triangles[i] = new int***[length_];
            for (int j = 0; j<length_; j++){
                //triangles[i][j] = new int[8][3][2];
                triangles[i][j] = Lattice::trianglesCCW_(i, j);
            }
        }
        triangles_ = triangles;
    }
    
    int*** Lattice::trianglesCCW_(int i, int j){
        int *** triangles = new int**[8];
        for (int n = 0; n < 8; n++){
            triangles[n] = new int*[3];
            for (int k = 0; k < 3; k++){
                triangles[n][k] = new int[2];
            }
        }
        
        //triangle 1 
        triangles[0][0][0] = i;
        triangles[0][0][1] = j;
        triangles[0][1][0] = Lattice::plusOne_(i);
        triangles[0][1][1] = Lattice::minusOne_(j);
        triangles[0][2][0] = Lattice::plusOne_(i);
        triangles[0][2][1] = j;

        //triangle 2
        triangles[1][0][0] = i;
        triangles[1][0][1] = j;
        triangles[1][1][0] = i;
        triangles[1][1][1] = Lattice::minusOne_(j);
        triangles[1][2][0] = Lattice::plusOne_(i);
        triangles[1][2][1] = Lattice::minusOne_(j);

        //triangle 3
        triangles[2][0][0] = i;
        triangles[2][0][1] = j;
        triangles[2][1][0] = Lattice::minusOne_(i);
        triangles[2][1][1] = Lattice::minusOne_(j);
        triangles[2][2][0] = i;
        triangles[2][2][1] = Lattice::minusOne_(j);

        //triangle 4
        triangles[3][0][0] = i;
        triangles[3][0][1] = j;
        triangles[3][1][0] = Lattice::minusOne_(i);
        triangles[3][1][1] = j;
        triangles[3][2][0] = Lattice::minusOne_(i);
        triangles[3][2][1] = Lattice::minusOne_(j);

        //triangle 5
        triangles[4][0][0] = i;
        triangles[4][0][1] = j;
        triangles[4][1][0] = Lattice::minusOne_(i);
        triangles[4][1][1] = Lattice::plusOne_(j);
        triangles[4][2][0] = Lattice::minusOne_(i);
        triangles[4][2][1] = j;

        //triangle 6
        triangles[5][0][0] = i;
        triangles[5][0][1] = j;
        triangles[5][1][0] = i;
        triangles[5][1][1] = Lattice::plusOne_(j);
        triangles[5][2][0] = Lattice::minusOne_(i);
        triangles[5][2][1] = Lattice::plusOne_(j);

        //triangle 7
        triangles[6][0][0] = i;
        triangles[6][0][1] = j;
        triangles[6][1][0] = Lattice::plusOne_(i);
        triangles[6][1][1] = Lattice::plusOne_(j);
        triangles[6][2][0] = i;
        triangles[6][2][1] = Lattice::plusOne_(j);

        //triangle 8
        triangles[7][0][0] = i;
        triangles[7][0][1] = j;
        triangles[7][1][0] = Lattice::plusOne_(i);
        triangles[7][1][1] = j;
        triangles[7][2][0] = Lattice::plusOne_(i);
        triangles[7][2][1] = Lattice::plusOne_(j);
        
        return triangles;
    }
    
    double Lattice::locQL_(int i, int j, int n, bool use_arccos){
        //Calculates QL on the nth triangle with central vertex i,j
        double rho, rho2, QLcos, QLsin;
        int i1 = triangles_[i][j][n][0][0];
        int j1 = triangles_[i][j][n][0][1];
        int i2 = triangles_[i][j][n][1][0];
        int j2 = triangles_[i][j][n][1][1];
        int i3 = triangles_[i][j][n][2][0];
        int j3 = triangles_[i][j][n][2][1];
        double *phi1 = Lattice::getPhi(i1,j1);
        double *phi2 = Lattice::getPhi(i2,j2);
        double *phi3 = Lattice::getPhi(i3,j3);
        rho2 = 2.*(1. + dot(phi1, phi2))*(1. + dot(phi2, phi3))*(1. + dot(phi3, phi1));
        rho = std::sqrt(rho2);
        QLcos = (1. + dot(phi1, phi2) + dot(phi2, phi3) + dot(phi3, phi1))/rho;
        QLsin = dot(phi1,cross(phi2,phi3))/rho;
#ifdef EXTREME_TESTING_MODE
        std:: cout << "QL from sine: "<< std::asin(QLsin)/(2.*M_PI) << std::endl;
        std:: cout << "QL from cosine: "<< std::acos(QLcos)/(2.*M_PI) << std::endl;
#endif
        if (use_arccos){ 
            return std::acos(QLcos)/(2.*M_PI);
        }
        else{
            return std::asin(QLsin)/(2.*M_PI);
        }
    }
    
    int* Lattice::getNeighbors_(int i, int j){
        static int nn[8];
        nn[0] = Lattice::plusOne_(i);
        nn[1] = j;
        nn[2] = i;
        nn[3] = Lattice::plusOne_(j);
        nn[4] = Lattice::minusOne_(i);
        nn[5] = j;
        nn[6] = i;
        nn[7] = Lattice::minusOne_(j);
        return nn;
    }
    
    double** Lattice::getNeighborPhis_(int i, int j){
        double ** nnPhis = new double*[4];
        nnPhis[0] = Lattice::getPhi(Lattice::plusOne_(i), j);
        nnPhis[1] = Lattice::getPhi(i, Lattice::plusOne_(j));
        nnPhis[2] = Lattice::getPhi(Lattice::minusOne_(i), j);
        nnPhis[3] = Lattice::getPhi(i, Lattice::minusOne_(j));
        return nnPhis;
    }
    
    void Lattice::printPhi_(int i, int j){
        double *phi = Lattice::getPhi(i,j);
        std::cout << "At point (" << i << "," << j <<"), phi = (";
        std::cout << phi[0] << "," << phi[1] << ","<< phi[2] << ")" << std::endl;
    }
    
}//end class definition
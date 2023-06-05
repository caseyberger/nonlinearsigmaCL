// Casey Berger
// Created: Mar 28, 2023
// Last edited: June 1, 2023

#include <iostream> //cout, endl
#include <cmath> //sqrt, sin, cos, acos, asin, exp, abs, remainder
#include "mathlib.h" //dot, cross

#include "lattice.h"


namespace nonlinearsigma{
    //public functions
    //constructor
    Lattice::Lattice(int length, double beta, double itheta){
        Lattice::setLength(length); //set length
        Lattice::setBeta(beta); //set beta
        Lattice::setiTheta(itheta); //set itheta
        fixedr_ = false; //leave rng untouched (this should only be true when testing)
    }
    //other public functions
    void Lattice::setLength(int length){
        //tested 6/1/2023
        //if you do this after initializing, you must reinitialize!
        length_ = length;
    }
    
    void Lattice::setBeta(double beta){
        //tested 6/1/2023
        beta_ = beta;
    }
    
    void Lattice::setiTheta(double itheta){
        //tested 6/1/2023
        itheta_ = itheta;
    }
    
    void Lattice::fixRNG(double r1, double r2){
        //tested 6/1/2023
        fixedr_ = true;
        r1_ = r1;
        r2_ = r2;
    }
    
    void Lattice::freeRNG(){
        //tested 6/1/2023
        fixedr_ = false;
    }
    
    int Lattice::getLength(){
        //tested 6/1/2023
        return length_;
    }
    
    double Lattice::getBeta(){
        //tested 6/1/2023
        return beta_;
    }
    
    double Lattice::getiTheta(){
        //tested 6/1/2023
        return itheta_;
    }
    
    double* Lattice::getPhi(int i, int j){
        //tested 5/30/2023
        return grid_[i][j];
    }
    
    double* Lattice::getRandNums(){
        //tested 6/1/2023
        static double r[2] = {r1_, r2_};
        return r;
    }
    
    double Lattice::getPhiMag(int i, int j){
        //tested 6/1/2023
        double *phi = Lattice::getPhi(i,j);
        double phi_mag = dot(phi,phi);
        return phi_mag;
    }
    
    double Lattice::getPhiTot(){
        //tested 6/1/2023
        double phi_tot = 0.;
        for(int i = 0; i < length_; i++){
            for(int j = 0; j < length_; j++){
                phi_tot += Lattice::getPhiMag(i, j);
            }
        }
        return phi_tot;
    }
    
    void Lattice::initialize(){
        //tested 5/30/2023
        double *** grid = new double**[length_];
        for(int i = 0; i < length_; i++){
            grid[i] = new double*[length_];
        }
        //allocation - 3 phi components
        for(int i = 0; i < length_; i++){
            for (int j = 0; j<length_; j++){
                grid[i][j] = new double[3];
                double *phi = Lattice::makePhi_();
                grid[i][j][0] = phi[0];
                grid[i][j][1] = phi[1];
                grid[i][j][2] = phi[2];
            }
        }
        grid_ = grid;
        Lattice::makeTriangles_();
        Lattice::zeroCount();
    }
    
    void Lattice::printLattice(){
        //tested 5/30/2023
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                Lattice::printPhi_(i, j);
            }
        }
    }
    
    void Lattice::printTriangles(int i, int j){
        //tested 6/1/2023
        std::cout << "At point (" << i << "," << j << ")," << std::endl;
        for (int n = 0; n < 8; n++){
            std::cout << "Triangle "<< n + 1 << " = (";
            std::cout << triangles_[i][j][n][0][0] << ","<< triangles_[i][j][n][0][1] << "), (";
            std::cout << triangles_[i][j][n][1][0] << ","<< triangles_[i][j][n][1][1] << "), (";
            std::cout << triangles_[i][j][n][2][0] << ","<< triangles_[i][j][n][2][1]<< ")" << std::endl;
        }
    }
    
    void Lattice::printNeighbors(int i, int j){
        //tested 6/1/2023
        int *nn = Lattice::getNeighbors_(i,j);
        double **nnphi = Lattice::getNeighborPhis_(i,j);
        std::cout << "At (i,j) = " << i << "," << j << " the neighbors are: " << std::endl;
        std::cout << "(" << nn[0] << "," << nn[1] << "), with phi (" << nnphi[0][0] << "," << nnphi[0][1] << "," << nnphi[0][2] << ")" << std::endl;
        std::cout << "(" << nn[2] << "," << nn[3] << "), with phi (" << nnphi[1][0] << "," << nnphi[1][1] << "," << nnphi[1][2] << ")" << std::endl;
        std::cout << "(" << nn[4] << "," << nn[5] << "), with phi (" << nnphi[2][0] << "," << nnphi[2][1] << "," << nnphi[2][2] << ")" << std::endl;
        std::cout << "(" << nn[6] << "," << nn[7] << "), with phi (" << nnphi[3][0] << "," << nnphi[3][1] << "," << nnphi[3][2] << ")" << std::endl;
    }
    
    double Lattice::calcQL(){
        //tested 6/1/2023
        //calculates topological charge -- note this does not produce integer values!
        //This Q_L is not renormalized. You can renormalize it later with Z 
        //is renormalizing what will make it an integer?
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
        return Q_L;
    }
    
    double Lattice::calcAL(){
        //tested 6/1/2023
        //calculates the standard lattice action A_L
        //you may be double counting things or you may be half counting. 
        //If you are off by 1/2 or 2, check here first
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
        //tested 6/1/2023
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
#ifdef EXTREME_TESTING_MODE
                Lattice::printPhi_(i, j);
                std:: cout << "Phi old = (" << phi_old[0] << "," << phi_old[1] << "," phi_old[2] << ")"  std::endl;
#endif
                //update lattice
                double *phi_new = Lattice::makePhi_();
                grid_[i][j] = phi_new;
                Lattice::printPhi_(i, j);
                std:: cout << "Phi old = (" << phi_new[0] << "," << phi_new[1] << "," phi_new[2] << ")"  std::endl;
                Sf = Lattice::calcSL();
                dS = Sf - Si;
#ifdef TEST_CONSTANT_RN
                r = 0.5;
#else
                r = ((double)std::rand())/((double)RAND_MAX);
#endif
                if(dS < 0 || std::exp(dS) > r){
                    acceptCount_++;//increment accept counter
#ifdef EXTREME_TESTING_MODE
                    std:: cout << "Accept" << std::endl;
#endif
                }
                else{
                    grid_[i][j] = phi_old;//change the value back to the old phi
                    rejectCount_++;//increment reject counter
#ifdef EXTREME_TESTING_MODE
                    std:: cout << "Reject" << std::endl;
#endif
                }
            }//loop over j
        }//loop over i
        double acc_rate = (double)acceptCount_/((double)acceptCount_ + (double)rejectCount_);
        accRate_ = acc_rate;
#ifdef EXTREME_TESTING_MODE
        std:: cout << "Acceptance rate: " << acc_rate << std::endl;
#endif
    }
    
    void Lattice::thermalize(int ntherm){

        for (int n = 0; n < ntherm; n++){
#ifdef EXTREME_TESTING_MODE
            std::cout << "Thermalization step " << n << std::endl;
#endif
            Lattice::metropolisStep();
        }
        //should you zero the count at the end of this?
    }
    
    void Lattice::zeroCount(){
        acceptCount_ = 0;
        rejectCount_ = 0;
        accRate_ = 0.;
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
        if(fixedr_){
            r1 = r1_;
            r2 = r2_;
        }
        else{
            r1 = ((double)std::rand())/((double)RAND_MAX);
            r2 = ((double)std::rand())/((double)RAND_MAX);
        }
#endif
        //generate a random polar and azimuthal angle
        double inclination =   std::acos(1. - 2. * r1); //polar angle = inclination
        //source: http://corysimon.github.io/articles/uniformdistn-on-sphere/
        double azimuth =  2. * M_PI * r2; //azimuthal angle = azimuth
        //create unit spin vector components from angles
        phi[0] = std::sin(inclination) * std::cos(azimuth);
        phi[1] = std::sin(inclination) * std::sin(azimuth);
        phi[2] = std::cos(inclination);
        
        r1_ = r1;//save for debugging
        r2_ = r2;//save for debugging
        return phi;
    }
    
    int Lattice::plusOne_(int i){
        //tested 5/30/2023
        if(i+1 == length_){ return 0;}
        else{return i+1;}
    }
    
    int Lattice::minusOne_(int i){
        //tested 5/30/2023
        if(i == 0){ return length_-1;}
        else{return i-1;}
    }
    
    void Lattice::makeTriangles_(){
        //tested 6/1/2023
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
        //tested 6/1/2023
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
        //tested 6/1/2023
        //Calculates QL on the nth triangle with central vertex i,j
        double rho, rho2, QLc, QLs;
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
        QLc = (1. + dot(phi1, phi2) + dot(phi2, phi3) + dot(phi3, phi1))/rho;
        QLs = dot(phi1,cross(phi2,phi3))/rho;
        double QLcos = std::acos(QLc)/(2.*M_PI);
        double QLsin = std::asin(QLs)/(2.*M_PI);
        if (use_arccos){ 
            //adjust arccos so it has the same domain as arcsin (-pi/2,pi/2)
            if (QLcos > 0.5*M_PI){
                return QLcos - 2*M_PI;
            }
            else if (QLcos = - QLsin){
                return -QLcos;
            }
            else{
                return QLcos;
            }
        }
        else{
            return QLsin;
        }
    }
    
    void Lattice::checkQL(int i, int j){
        //tested 6/1/2023
        //ensures we get the same QLtri with either cosine or sine
        //also ensures that QLtri is in the correct range of [-pi/2, pi/2]
        //also ensures that QL over all triangles is an integer (w/in e-15)
        double QLtot = 0;
        double tol = 1.0e-15;
        for (int n = 0; n < 8; n++){
            double QLcos, QLsin;
            bool use_cosine = true;
            bool use_sine = false;
            QLcos = Lattice::locQL_(i,j,n,use_cosine);
            QLsin = Lattice::locQL_(i,j,n,use_sine);
            if (QLcos != QLsin){
                std::cout << "QLcos = " << QLcos << ", QLsin = " << QLsin << std::endl;
            }
            if (QLsin > 0.5*M_PI || QLsin < -0.5*M_PI){
                std::cout << "QL of triangle outside range [-pi/2,pi/2]: " << QLsin << std::endl;
            }
            QLtot += QLsin;
        }
        std::cout << "QLtot = " << QLtot << std::endl;
        if (std::abs(std::remainder(QLtot,1)) > tol){
            std::cout << "QL not an integer value: " << QLtot << std::endl;
        }
    }
    
    int* Lattice::getNeighbors_(int i, int j){
        //tested 6/1/2023
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
        //tested 6/1/2023
        double ** nnPhis = new double*[4];
        nnPhis[0] = Lattice::getPhi(Lattice::plusOne_(i), j);
        nnPhis[1] = Lattice::getPhi(i, Lattice::plusOne_(j));
        nnPhis[2] = Lattice::getPhi(Lattice::minusOne_(i), j);
        nnPhis[3] = Lattice::getPhi(i, Lattice::minusOne_(j));
        return nnPhis;
    }
    
    void Lattice::printPhi_(int i, int j){
        //tested 5/30/2023
        double *phi = Lattice::getPhi(i,j);
        std::cout << "At point (" << i << "," << j <<"), phi = (";
        std::cout << phi[0] << "," << phi[1] << ","<< phi[2] << ")" << std::endl;
    }
    
}//end class definition
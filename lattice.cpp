// Casey Berger
// Created: Mar 28, 2023
// Last edited: June 16, 2023

#include <iostream> //cout, endl
#include <cmath> //sqrt, sin, cos, acos, asin, exp, abs, remainder
#include <string> //string
#include <vector> //vector
#include <numeric> // iota
#include <algorithm>  // shuffle
#include <random> //default_random_engine
#include <array> 
#include "mathlib.h" //dot, cross
#include "lattice.h"


namespace nonlinearsigma{
    //public functions
    //constructor
    Lattice::Lattice(int length, double beta, double itheta){
        Lattice::setLength(length); //set length
        Lattice::setBeta(beta); //set beta
        Lattice::setiTheta(itheta); //set itheta
        Lattice::setnTherm(1000); //set therm steps to default number
        Lattice::setnMC(1000); //set Monte Carlo steps to default number
        
        Lattice::generateFilename_();
        fixedr_ = false; //this should only be set to true when testing
    }
    //other public functions
    void Lattice::setLength(int length){
        //tested 6/1/2023
        //if you do this after initializing, you must reinitialize!
        length_ = length;
        Lattice::generateFilename_();
    }
    
    void Lattice::setBeta(double beta){
        //tested 6/1/2023
        beta_ = beta;
        Lattice::generateFilename_();
    }
    
    void Lattice::setiTheta(double itheta){
        //tested 6/1/2023
        itheta_ = itheta;
        Lattice::generateFilename_();
    }
    
    void Lattice::setPhi(int i, int j, Lattice::field phi){
        //tested 6/5/2023
        grid_[i][j][0] = phi[0];
        grid_[i][j][1] = phi[1];
        grid_[i][j][2] = phi[2];
    }
    
    void Lattice::setAvgG(int i, int j, double Gij){
        Gij_[i][j]= Gij;
    }
    
    void Lattice::setnTherm(int ntherm){
        nTherm_ = ntherm;
        Lattice::generateFilename_();
    }
    
    void Lattice::setnMC(int nMC){
        nMC_ = nMC;
        Lattice::generateFilename_();
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
    
    int Lattice::getnTherm(){
       return nTherm_; 
    }
    
    int Lattice::getnMC(){
        return nMC_;
    }
    
    std::string Lattice::getFilename(){
        return filename_;
    }
    
    field Lattice::getPhi(int i, int j){
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
        field phi = Lattice::getPhi(i,j);
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
    
    double Lattice::getAvgG(int i, int j){
        double Gij = Gij_[i][j];
        return Gij;
    }
    
    void Lattice::initialize(){
        //tested 5/30/2023
        std::vector < std::vector < Lattice::field > > grid;
        std::vector < std::vector < double > > Gij;
        for(int i = 0; i < length_; i++){
            std::vector <double> Gj;
            std::vector < Lattice::field > gridj;
            for (int j = 0; j<length_; j++){
                field phi = Lattice::makePhi_();
                gridj.push_back(phi);
                Gj.push_back(0.);
            }
            Gij.push_back(Gj);
            grid.push_back(gridj);
        }
        grid_ = grid;
        Gij_ = Gij;
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
        std::array < Lattice::vertex, 4 > nn = Lattice::getNeighbors_(i,j);
        std::array < Lattice::field, 4 > nnphi = Lattice::getNeighborPhis_(i,j);
        std::cout << "At (i,j) = " << i << "," << j << " the neighbors are: " << std::endl;
        std::cout << "(" << nn[0][0] << "," << nn[0][1] << "), with phi (" << nnphi[0][0] << "," << nnphi[0][1] << "," << nnphi[0][2] << ")" << std::endl;
        std::cout << "(" << nn[1][0] << "," << nn[1][1] << "), with phi (" << nnphi[1][0] << "," << nnphi[1][1] << "," << nnphi[1][2] << ")" << std::endl;
        std::cout << "(" << nn[2][0] << "," << nn[2][1] << "), with phi (" << nnphi[2][0] << "," << nnphi[2][1] << "," << nnphi[2][2] << ")" << std::endl;
        std::cout << "(" << nn[3][0] << "," << nn[3][1] << "), with phi (" << nnphi[3][0] << "," << nnphi[3][1] << "," << nnphi[3][2] << ")" << std::endl;
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
                Lattice::field phi = Lattice::getPhi(i,j);
                std::array < Lattice::field, 4> phiNN = Lattice::getNeighborPhis_(i,j); //0 and 1 are + direction
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
    
    double Lattice::twoPointG(int i, int j){
        Lattice::field phi_00 = Lattice::getPhi(0,0);
        Lattice::field phi_ij = Lattice::getPhi(i,j);
        double G = 0;
        G = dot(phi_00, phi_ij);
        double oldAvgG = Lattice::getAvgG(i, j);
        int n = acceptCount_+rejectCount_;
        double newAvgG = (oldAvgG*n)/(n+1)+G/(n+1);
        Lattice::setAvgG(i,j,newAvgG);
        return G;
    }
    
    double Lattice::calcXi(){
        double Xi = 0.;
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                Xi += twoPointG(i, j);
            }
        }
        return Xi;
    }
    
    double* Lattice::calcF(){
        double F_Re = 0.;
        double F_Im = 0.;
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                F_Re += 0.5*twoPointG(i, j)*(std::cos(2.*M_PI*i/length_) + std::cos(2.*M_PI*j/length_));
                F_Im += 0.5*twoPointG(i, j)*(std::sin(2.*M_PI*i/length_) + std::sin(2.*M_PI*j/length_));
            }
        }
        static double F[2] = {F_Re, F_Im};
        return F;
    }
    
    
    void Lattice::metropolisStep(){
        //tested 6/5/2023
        double Si, Sf, dS, r;
        
        int nsites = length_*length_;
        std::vector<int> site_arr(nsites);
        std::iota(site_arr.begin(), site_arr.end(), 0);     
        shuffle(site_arr.begin(), site_arr.end(), std::default_random_engine(1232));

        for(unsigned int n = 0; n < site_arr.size(); n++){
            int i = site_arr[n]/length_;
            int j = site_arr[n]%length_;
            Si = Lattice::calcSL();
            Lattice::field phi_old = Lattice::getPhi(i, j);

            //update lattice
            Lattice::field phi_new = Lattice::makePhi_();
            Lattice::setPhi(i, j, phi_new);
            Sf = Lattice::calcSL();
            dS = Sf - Si;
#ifdef TEST_CONSTANT_RN
            r = 0.5;
#else
            r = ((double)std::rand())/((double)RAND_MAX);
#endif
            if(dS < 0 || r < std::exp(-1.*dS)){
                acceptCount_++;//increment accept counter
            }
            else{
                Lattice::setPhi(i, j, phi_old);//change the value back to the old phi
                rejectCount_++;//increment reject counter
            }
        }//loop over sites
        double acc_rate = (double)acceptCount_/((double)acceptCount_ + (double)rejectCount_);
        accRate_ = acc_rate;
#ifdef EXTREME_TESTING_MODE
        std:: cout << "Acceptance rate: " << acc_rate << std::endl;
#endif
        site_arr.clear();
    }
    
    void Lattice::thermalize(){
        //tested 6/5/2023
        for (int n = 0; n < nTherm_; n++){
#ifdef EXTREME_TESTING_MODE
            std::cout << "Thermalization step " << n << std::endl;
#endif
            Lattice::metropolisStep();
        }
        Lattice::zeroCount();
    }
    
    void Lattice::zeroCount(){
        //tested 6/5/2023
        acceptCount_ = 0;
        rejectCount_ = 0;
        accRate_ = 0.;
    }
    
    double Lattice::acceptanceRate(){
        //tested 6/5/2023
        return accRate_;
    }
    
    //private functions
    Lattice::field Lattice::makePhi_(){
        //tested 6/1/2023
        Lattice::field phi;
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
        //tested 6/16/2023
        
        std::vector < std::vector < site_triangles > > all_triangles;
        for(int i = 0; i < length_; i++){
            std::vector < site_triangles > tri_y;
            for (int j = 0; j<length_; j++){
                site_triangles local_triangles;
                //triangle 1
                Lattice::vertex v1 = {i,j};
                Lattice::vertex v2 = {Lattice::plusOne_(i),Lattice::minusOne_(j)};
                Lattice::vertex v3 = {Lattice::plusOne_(i),j};
                local_triangles[0] = {v1, v2, v3};
                
                //triangle 2
                v3 = v2;
                v2 = {i,Lattice::minusOne_(j)};
                local_triangles[1] = {v1, v2, v3};

                //triangle 3
                v3 = v2;
                v2 = {Lattice::minusOne_(i),Lattice::minusOne_(j)};
                local_triangles[2] = {v1, v2, v3};
                              
                //triangle 4
                v3 = v2;
                v2 = {Lattice::minusOne_(i),j};
                local_triangles[3] = {v1, v2, v3};

                //triangle 5
                v3 = v2;
                v2 = {Lattice::minusOne_(i),Lattice::plusOne_(j)};
                local_triangles[4] = {v1, v2, v3};
                
                //triangle 6
                v3 = v2;
                v2 = {i,Lattice::plusOne_(j)};
                local_triangles[5] = {v1, v2, v3};

                //triangle 7
                v3 = v2;
                v2 = {Lattice::plusOne_(i),Lattice::plusOne_(j)};
                local_triangles[6] = {v1, v2, v3};

                //triangle 8
                v3 = v2; 
                v2 = {Lattice::plusOne_(i),j};
                local_triangles[7] = {v1, v2, v3};
                
                tri_y.push_back(local_triangles);
            }
            all_triangles.push_back(tri_y);
        }
        triangles_ = all_triangles;
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
        Lattice::field phi1 = Lattice::getPhi(i1,j1);
        Lattice::field phi2 = Lattice::getPhi(i2,j2);
        Lattice::field phi3 = Lattice::getPhi(i3,j3);
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
    
    std::array < Lattice::vertex, 4> Lattice::getNeighbors_(int i, int j){
        //tested 6/1/2023
        std::array < Lattice::vertex, 4> nn;
        nn[0][0] = Lattice::plusOne_(i);
        nn[0][1] = j;
        nn[1][0] = i;
        nn[1][1] = Lattice::plusOne_(j);
        nn[2][0] = Lattice::minusOne_(i);
        nn[2][1] = j;
        nn[3][0] = i;
        nn[3][1] = Lattice::minusOne_(j);
        return nn;
    }
    
    std::array < Lattice::field, 4 > Lattice::getNeighborPhis_(int i, int j){
        //tested 6/1/2023
        std::array < Lattice::field, 4> nnPhis;
        nnPhis[0] = Lattice::getPhi(Lattice::plusOne_(i), j);
        nnPhis[1] = Lattice::getPhi(i, Lattice::plusOne_(j));
        nnPhis[2] = Lattice::getPhi(Lattice::minusOne_(i), j);
        nnPhis[3] = Lattice::getPhi(i, Lattice::minusOne_(j));
        return nnPhis;
    }
    
    void Lattice::printPhi_(int i, int j){
        //tested 5/30/2023
        Lattice::field phi = Lattice::getPhi(i,j);
        std::cout << "At point (" << i << "," << j <<"), phi = (";
        std::cout << phi[0] << "," << phi[1] << ","<< phi[2] << ")" << std::endl;
    }
    
    void Lattice::generateFilename_(){
        std::string l_str   = std::to_string(length_);
        std::string b_str   = std::to_string(beta_);
        std::string th_str  = std::to_string(itheta_);
        std::string nt_str  = std::to_string(nTherm_);
        std::string nmc_str = std::to_string(nMC_);
        std::string fname = "nonlinearsigma_data_L_" + l_str + "_beta_" + b_str + "_itheta_" + th_str + "_ntherm_" + nt_str + "_nMC_" + nmc_str + ".csv";
        filename_ = fname;
    }
    
}//end class definition
// Casey Berger
// Created: Mar 28, 2023
// Last edited: Feb 16, 2023
/* changelog for last edit: 
 - added new function exceptionalTriangles to check all six triangles that touch vertex (i,j)
 - updated removeExceptional to use the new exceptionalTriangles function
 
 
 suggestions for changes
 - can you create a class in C++ that inherits the lattice object? 
         If so, make a testing object so you can put functions like fixRNG 
        and freeRNG and setTrig into another class to streamline this one.
*/
#include <iostream> //cout, endl
#include <cmath> //sqrt, sin, cos, acos, asin, exp, abs, remainder
#include <string> //string
#include <vector> //vector
#include <numeric> // iota
#include <algorithm>  // shuffle
#include <random> //default_random_engine
#include <array> 
#include <omp.h>
#include <fstream> //fout

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
        Lattice::setFreq(100); //set frequency between saved configs to default number
        Lattice::generateFilename_();
        fixedr_ = false; //this should only be set to true when testing
        use_arcsin_ = true;
    }
    //other public functions
    void Lattice::setLength(int length){
        //tested 6/1/2023
        length_ = length;
        Lattice::generateFilename_();
        Lattice::initialize();
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
    
    void Lattice::setnTherm(int ntherm){
        nTherm_ = ntherm;
        Lattice::generateFilename_();
    }
    
    void Lattice::setnMC(int nMC){
        nMC_ = nMC;
        Lattice::generateFilename_();
    }
    
    void Lattice::setFreq(int freq){
        freq_ = freq;
        Lattice::generateFilename_();
    }
    
    void Lattice::setTrig(bool use_arcsin){
        use_arcsin_ = use_arcsin;
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
    
    Lattice::field Lattice::getPhi(int i, int j){
        //tested 5/30/2023
        return grid_[i][j];
    }
    
    int Lattice::getAttempts(int i, int j){
        return gridAttempts_[i][j];
    }
    
    double* Lattice::getRandNums(){
        //tested 6/1/2023
        static double r[2] = {r1_, r2_};
        return r;
    }
    
    double Lattice::getPhiTot(){
        //tested 6/1/2023
        //double phi_tot = 0.;
        double phi_tot(0.);//optimization 7/4/23
        #pragma omp parallel for collapse(2) default(none) shared(length_) reduction(+:phi_tot)
        for(int i = 0; i < length_; i++){
            for(int j = 0; j < length_; j++){
                field phi(Lattice::getPhi(i,j)); 
                phi_tot += dot(phi,phi);
            }
        }
        return phi_tot;
    }
    
    double Lattice::getAvgG(int i, int j){
        //double Gij = Gij_[i][j];
        double Gij(Gij_[i][j]);//optimization 7/4/23
        return Gij;
    }
    
    void Lattice::removeExceptional(int i, int j, int exc_lim){
        bool exceptional_config = true;
        int exc_count = 0;
        while (exceptional_config){
            //checks the six triangles that inclue vertex (i,j) 
            if (Lattice::exceptionalTriangles(i, j, exc_count)){
                exceptional_config = true;
                exc_count++;
                //try new field value
                field phi_new = Lattice::makePhi_();
                grid_[i][j][0] = phi_new[0];
                grid_[i][j][1] = phi_new[1];
                grid_[i][j][2] = phi_new[2]; 
            }
            else{
                exceptional_config = false;
            }
            if (exc_count >= exc_lim){
                break;
            }
        }//while still exceptional at vertex (i,j)
        gridAttempts_[i][j] = exc_count;//update grid with number of attempts made
    }
    
    bool Lattice::exceptionalTriangles(int i, int j, int exc_count){
        //triangles owned by vertex i,j
        //algorithm: check each triangle individually and print (in test mode) which is exceptional
        bool any_exceptional = false;
        bool tri1, tri2, tri3, tri4, tri5, tri6;
        int im1 = Lattice::minusOne_(i);
        int jm1 = Lattice::minusOne_(j);
        //triangle 1: vertex (i,j), lower triangle
        if (Lattice::exceptionalConfig(i,j,0)){
            tri1 = true;
            any_exceptional = true;
            }
        else{tri1 = false;}
        //triangle 2: vertex (i,j), upper triangle
        if (Lattice::exceptionalConfig(i,j,1)){
            tri2 = true;
            any_exceptional = true;
        }
        else{tri2 = false;}
        //triangle 3: vertex (i-1,j-1), lower triangle
        if (Lattice::exceptionalConfig(im1,jm1,0)){
            tri3 = true;
            any_exceptional = true;
        }
        else{tri3 = false;}
        //triangle 4: vertex (i-1,j-1), upper triangle
        if (Lattice::exceptionalConfig(im1,jm1,1)){
            tri4 = true;
            any_exceptional = true;
        }
        else{tri4 = false;}
        //triangle 5: vertex (i-1,j), lower triangle
        if (Lattice::exceptionalConfig(im1,j,0)){
            tri5 = true;
            any_exceptional = true;
        }
        else{tri5 = false;}
        //triangle 6: vertex (i,j-1), upper triangle
        if (Lattice::exceptionalConfig(i,jm1,1)){
            tri6 = true;
            any_exceptional = true;
        }
        else{tri6 = false;}
        
        
#ifdef TESTING_MODE
        if (any_exceptional){
            Lattice::printPhi_(i, j);
            std::cout << "Attempt "<< exc_count << ":" << std::endl;
            if (tri1){
                std::cout << "Triangle 1 exceptional "<<std::endl; 
            }
            if (tri2){
                std::cout << "Triangle 2 exceptional "<<std::endl;
            }
            if (tri3){
                std::cout << "Triangle 3 exceptional "<<std::endl;
            }
            if (tri4){
                std::cout << "Triangle 4 exceptional "<<std::endl;
            }
            if (tri5){
                std::cout << "Triangle 5 exceptional "<<std::endl;
            }
            if (tri6){
                std::cout << "Triangle 6 exceptional "<<std::endl;
            }
        }
#endif
        return any_exceptional;
    }
    
    
    void Lattice::clean(){
        //cleaning the lattice means removing exceptional configurations
        int nsites(length_*length_); 
        std::vector<int> site_arr(nsites);
        std::iota(site_arr.begin(), site_arr.end(), 0);     
        shuffle(site_arr.begin(), site_arr.end(), std::default_random_engine(1232));
        
        int exc_lim(10000);

        for(unsigned int n = 0; n < site_arr.size(); n++){
            int i(site_arr[n]/length_);
            int j(site_arr[n]%length_);
            Lattice::removeExceptional(i, j, exc_lim);
        }//loop over sites
    }
    
    void Lattice::initialize(){
        //tested 5/30/2023
        std::vector < std::vector < Lattice::field > > grid;
        std::vector < std::vector < double > > Gij;
        std::vector < std::vector < int > > gridAttempts;
        std::vector < std::vector < bool > > gridMCAccepted;
        for(int i = 0; i < length_; i++){
            std::vector <double> Gj;
            std::vector <int> gridAttemptsj;
            std::vector <bool> gridMCAcceptedj;
            std::vector < Lattice::field > gridj;
            for (int j = 0; j<length_; j++){
                field phi = Lattice::makePhi_();
                gridj.push_back(phi);
                Gj.push_back(0.);
                gridAttemptsj.push_back(0);
                gridMCAcceptedj.push_back(false);
            }
            Gij.push_back(Gj);
            grid.push_back(gridj);
            gridAttempts.push_back(gridAttemptsj);
            gridMCAccepted.push_back(gridMCAcceptedj);
        }
        grid_ = grid;
        Gij_ = Gij;
        gridAttempts_ = gridAttempts;
        gridMCAccepted_ = gridMCAccepted;
        Lattice::makeTriangles_();
        Lattice::zeroCount();
        Lattice::clean();
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
        //updated 7/14/2023
        std::cout << "At point (" << i << "," << j << ")," << std::endl;
        std::cout << "Triangle 1 = (";
        std::cout << triangles_[i][j][0][0][0] << ","<< triangles_[i][j][0][0][1] << "), (";
        std::cout << triangles_[i][j][0][1][0] << ","<< triangles_[i][j][0][1][1] << "), (";
        std::cout << triangles_[i][j][0][2][0] << ","<< triangles_[i][j][0][2][1]<< ")" << std::endl;
        std::cout << "Triangle 2 = (";
        std::cout << triangles_[i][j][1][0][0] << ","<< triangles_[i][j][1][0][1] << "), (";
        std::cout << triangles_[i][j][1][1][0] << ","<< triangles_[i][j][1][1][1] << "), (";
        std::cout << triangles_[i][j][1][2][0] << ","<< triangles_[i][j][1][2][1]<< ")" << std::endl;
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
    
    void Lattice::saveConfig(int step){
        std::string step_str   = std::to_string(step);
        std::string fname = "config_"+step_str+".csv";
        std::ofstream fout; //output stream
        fout.open(fname.c_str(),std::ios::out);

        // check if files are open
        if (!fout.is_open())
        {
            std::cerr << "Unable to open file " << fname <<"." << std::endl;
            std::exit(10);
        }
        fout.setf(std::ios::fixed);
        fout << "i,j,phi_x,phi_y,phi_z,Gij, numAttempts, exceptional, MCAccepted" << std::endl;
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                field phi = Lattice::getPhi(i,j);
                double Gij = Lattice::getAvgG(i,j);
                int numAttempts = Lattice::getAttempts(i,j);
                fout << i <<","<< j << ",";
                fout << phi[0] << "," << phi[1] << "," << phi[2] << ",";
                fout << Gij << "," <<  numAttempts;
                if (Lattice::exceptionalConfig(i,j,0) or Lattice::exceptionalConfig(i,j,1)){fout << ", Y";}
                else{fout << ", N";}
                if (gridMCAccepted_[i][j]){fout << ", Y" << std::endl;}
                else{fout << ", N" << std::endl;}
            }
        }
        fout.close();
    }
    
    double Lattice::calcQL(){
        //updated 7/14/2023 to deal with new triangles
        //calculates topological charge
        double Q_L(0.);//optimization 7/4/23
        #pragma omp parallel for collapse(2) default(none) shared(length_,use_arcsin_) reduction(+:Q_L)
        for (int i = 0; i<length_; i++){
            for (int j = 0; j<length_; j++){
                //add triangle 1 at this vertex
                Q_L += Lattice::locQL_(i, j, 0, use_arcsin_);
                //add triangle 2 at this vertex
                Q_L += Lattice::locQL_(i, j, 1, use_arcsin_);
            }//loop over j
        }//loop over i
        return Q_L;
    }
    
    double Lattice::calcAL(){
        //tested 6/1/2023
        //calculates the standard lattice action A_L
        //you may be double counting things or you may be half counting. 
        //If you are off by 1/2 or 2, check here first
        double A_L(0.);
        #pragma omp parallel for collapse(2) default(none) shared(length_) reduction(+:A_L)
        for (int i = 0; i<length_; i++)
        {
            for (int j = 0; j<length_; j++)
            {
                Lattice::field phi(Lattice::getPhi(i,j));
                std::array < Lattice::field, 4> phiNN(Lattice::getNeighborPhis_(i,j));
                //nearest neighbors in positive direction:
                A_L += dot(phi, phiNN[0]) + dot(phi, phiNN[1]);
                //nearest neighbors in negative direction:
                //A_L += dot(phi, phiNN[2]) + dot(phi, phiNN[2]);
            }
        }
        return -1.*beta_*A_L;
    }
    
    double Lattice::calcSL(){
        //tested 6/1/2023
        //calculates the full lattice action S_L = A_L - i theta Q_L
        //optimization target -- remove this function
        double S_L = Lattice::calcAL() - 1.*itheta_*Lattice::calcQL();
        return S_L;
    }
    
    void Lattice::calcGij(){
        #pragma omp parallel for collapse(2) default(none) shared(length_)
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                Lattice::field phi_00(Lattice::getPhi(0,0));
                Lattice::field phi_ij(Lattice::getPhi(i,j));
                double G(dot(phi_00, phi_ij));
                double oldAvgG(Lattice::getAvgG(i, j));
                int n(acceptCount_+rejectCount_);
                double newAvgG = (oldAvgG*n)/(n+1)+G/(n+1);
                Gij_[i][j] = newAvgG;
            }
        }
    }
    
    double Lattice::calcXi(){
        //double Xi = 0.;
        double Xi(0.);//optimization 7/4/23
        #pragma omp parallel for collapse(2) default(none) shared(length_) reduction(+:Xi)
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                Xi += Lattice::getAvgG(i,j);
            }
        }
        return Xi;
    }
    
    double* Lattice::calcF(){
        //double F_Re = 0.;
        //double F_Im = 0.;
        double F_Re(0.),F_Im(0.);//optimization 7/4/23
        #pragma omp parallel for collapse(2) default(none) shared(length_) reduction(+:F_Re,F_Im)
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                double Gij = Lattice::getAvgG(i,j);
                F_Re += 0.5*Gij*(std::cos(2.*M_PI*i/length_) + std::cos(2.*M_PI*j/length_));
                F_Im += 0.5*Gij*(std::sin(2.*M_PI*i/length_) + std::sin(2.*M_PI*j/length_));
            }
        }
        static double F[2] = {F_Re, F_Im};
        return F;
    }
    
    
    void Lattice::metropolisStep(){
        //tested 6/5/2023
        double Si, Sf, dS, r;
        Lattice::field phi_old, phi_new;
        int exc_lim = 10000;
        
        //generate randomized array of sites
        int nsites(length_*length_); 
        std::vector<int> site_arr(nsites);
        std::iota(site_arr.begin(), site_arr.end(), 0);     
        shuffle(site_arr.begin(), site_arr.end(), std::default_random_engine(1232));

        for(unsigned int n = 0; n < site_arr.size(); n++){
            int i(site_arr[n]/length_);
            int j(site_arr[n]%length_);
            Si = Lattice::calcAL() - 1.*itheta_*Lattice::calcQL();
            phi_old = Lattice::getPhi(i, j);//store old field in case you reject the change
            
            //generate new field at site
            phi_new = Lattice::makePhi_();
            grid_[i][j][0] = phi_new[0];
            grid_[i][j][1] = phi_new[1];
            grid_[i][j][2] = phi_new[2];
            
            //make sure the change isn't exceptional
            Lattice::removeExceptional(i, j, exc_lim);
            
            
            //calculate change in action with new field 
            Sf = Lattice::calcAL() - 1.*itheta_*Lattice::calcQL();
            dS = Sf - Si;
#ifdef TEST_CONSTANT_RN
            r = 0.5;
#else
            r = ((double)std::rand())/((double)RAND_MAX);
#endif
            if(dS < 0 || r < std::exp(-1.*dS)){
                acceptCount_++;//increment accept counter
                gridMCAccepted_[i][j] = true; //update MC results for this sweep
            }
            else{
                grid_[i][j][0] = phi_old[0];
                grid_[i][j][1] = phi_old[1];
                grid_[i][j][2] = phi_old[2];
                rejectCount_++;//increment reject counter
                gridAttempts_[i][j] = 0; //reset num attempts to zero because we stuck with the old config
                gridMCAccepted_[i][j] = false; //update MC results for this sweep
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
    
    bool Lattice::exceptionalConfig(int i, int j, int n){
        //check that these are done counterconclockwise
        int i1(triangles_[i][j][n][0][0]);
        int j1(triangles_[i][j][n][0][1]);
        int i2(triangles_[i][j][n][1][0]);
        int j2(triangles_[i][j][n][1][1]);
        int i3(triangles_[i][j][n][2][0]);
        int j3(triangles_[i][j][n][2][1]);
        Lattice::field phi1(Lattice::getPhi(i1,j1));
        Lattice::field phi2(Lattice::getPhi(i2,j2));
        Lattice::field phi3(Lattice::getPhi(i3,j3));
        double check1 = dot(phi1,cross(phi2,phi3));
        double check2 = 1. + dot(phi1, phi2) + dot(phi2, phi3) + dot(phi3, phi1);
        if (check1 == 0 or check2 <= 0.){return true;}
        else {return false;}
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
        
        r1_ = r1;//save for debugging -- consider wrapping in a #ifndef statement for optimization
        r2_ = r2;//save for debugging -- consider wrapping in a #ifndef statement for optimization
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
        //updated 7/14/2023
        
        std::vector < std::vector < site_triangles > > all_triangles;
        for(int i = 0; i < length_; i++){
            std::vector < site_triangles > tri_y;
            for (int j = 0; j<length_; j++){
                site_triangles local_triangles;
                //triangle 1
                Lattice::vertex v1({i,j});
                Lattice::vertex v2({Lattice::plusOne_(i),Lattice::plusOne_(j)});
                Lattice::vertex v3({i,Lattice::plusOne_(j)});
                local_triangles[0] = {v1, v2, v3};
                
                //triangle 2
                v3 = v2;
                v2 = {Lattice::plusOne_(i),j};
                local_triangles[1] = {v1, v2, v3};
                
                tri_y.push_back(local_triangles);
            }
            all_triangles.push_back(tri_y);
        }
        triangles_ = all_triangles;
    }
    
    double Lattice::locQL_(int i, int j, int n, bool use_arcsin){
        //updated 7/14/2023 for new triangles
        //Calculates QL on the nth triangle with central vertex i,j
        double tol = 1.0e-5;
        double rho, QLc, QLs;
        int i1(triangles_[i][j][n][0][0]);
        int j1(triangles_[i][j][n][0][1]);
        int i2(triangles_[i][j][n][1][0]);
        int j2(triangles_[i][j][n][1][1]);
        int i3(triangles_[i][j][n][2][0]);
        int j3(triangles_[i][j][n][2][1]);
        Lattice::field phi1(Lattice::getPhi(i1,j1));
        Lattice::field phi2(Lattice::getPhi(i2,j2));
        Lattice::field phi3(Lattice::getPhi(i3,j3));
        rho = std::sqrt(2.*(1. + dot(phi1, phi2))*(1. + dot(phi2, phi3))*(1. + dot(phi3, phi1)));
        QLc = (1. + dot(phi1, phi2) + dot(phi2, phi3) + dot(phi3, phi1))/rho;
        QLs = dot(phi1,cross(phi2,phi3))/rho;
        double QLcos = std::acos(QLc)/(2.*M_PI);
        double QLsin = std::asin(QLs)/(2.*M_PI);
        if(use_arcsin){return QLsin;}
        else{//adjust arccos so it has the same domain as arcsin (-pi/2,pi/2)
            if (QLcos > 0.5*M_PI){QLcos += -2.*M_PI;}
            if (std::abs(QLcos + QLsin) < tol){QLcos *= -1.;}
            return QLcos;
        }
    }
    
    void Lattice::checkQL(){
        //tested 6/1/2023
        //ensures we get the same QLtri with either cosine or sine
        //also ensures that QLtri is in the correct range of [-pi/2, pi/2]
        //also ensures that QL over all triangles is an integer (w/in some tolerance)
        double QLtot = 0;
        double tol = 1.0e-5;
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                for (int n = 0; n < 2; n++){
                    double QLcos, QLsin;
                    bool use_cosine = true;
                    bool use_sine = false;
                    QLcos = Lattice::locQL_(i,j,n,use_cosine);
                    QLsin = Lattice::locQL_(i,j,n,use_sine);
                    //adjust arccos to match arcsin domain
                    if (QLcos > 0.5*M_PI){QLcos -= 2*M_PI;}
                    if (std::abs(QLcos + QLsin) < tol){QLcos = -1.*QLcos;}
                    //check if they are equivalent
                    if (std::abs(QLcos - QLsin) > tol){
                        std::cout << "QLcos = " << QLcos << ", QLsin = " << QLsin << std::endl;
                    }
                    else{
                        std::cout << "QLcos = QLsin = " << QLsin << std::endl;
                    }
                    if (QLsin > 0.5 || QLsin < -0.5){
                        std::cout << "QL of triangle outside range [-1/2,1/2]: " << QLsin << std::endl;
                    }
                    QLtot += QLsin;
                }//loop over n
            }//loop over j
        }//loop over i
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
        std::string frq_str = std::to_string(freq_);
        std::string fname = "nonlinearsigma_data_L_" + l_str + "_beta_" + b_str + "_itheta_" + th_str + "_ntherm_" + nt_str + "_nMC_" + nmc_str + "_freq_" + frq_str + ".csv";
        filename_ = fname;
    }
    
}//end class definition
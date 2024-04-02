// Casey Berger
// Created: Apr 2, 2024
// Last edited: Apr 2, 2024
/* 
Plans: this class will be the base class, and individual CL projects will inherit it
It should be able to do CL evolution if fed the drift function
It should store complex fields and numbers and do complex dot and cross products
*/
#include <iostream> //cout, endl
#include <cmath> //sqrt, sin, cos, acos, asin, exp, abs, remainder
#include <string> //string
#include <vector> //vector
#include <array> 
#include <omp.h>
#include <fstream> //fout
#include <numeric> // iota
#include <algorithm>  // shuffle
#include <random> //default_random_engine


#include "mathlib.h" //dot, cross
#include "cllib.h" //complex dot, complex cross (put those in here??)


namespace cl{
    //public functions
    //constructor
    ComplexLattice::Lattice(int length){
        ComplexLattice::SetLength(length); //set length
        ComplexLattice::GenerateFilename_();//set output filename
    }
    //other public functions
    void ComplexLattice::SetLength(int length){
        length_ = length;
        Lattice::generateFilename_();
        Lattice::initialize();
    }
    
    int ComplexLattice::GetLength(){
        //tested 6/1/2023
        return length_;
    }
    
    std::string ComplexLattice::GetFilename(){
        return filename_;
    }
    
    ComplexLattice::complex_field ComplexLattice::GetField(int i, int j){
        return lattice_[i][j];
    }
    
   
    void ComplexLattice::Initialize(){
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
    
    void ComplexLattice::PrintLattice(){
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                Lattice::printPhi_(i, j);
            }
        }
    }
    
    
    void ComplexLattice::PrintNeighbors(int i, int j){
        std::array < Lattice::vertex, 4 > nn = Lattice::getNeighbors_(i,j);
        std::array < Lattice::field, 4 > nnphi = Lattice::getNeighborPhis_(i,j);
        std::cout << "At (i,j) = " << i << "," << j << " the neighbors are: " << std::endl;
        std::cout << "(" << nn[0][0] << "," << nn[0][1] << "), with phi (" << nnphi[0][0] << "," << nnphi[0][1] << "," << nnphi[0][2] << ")" << std::endl;
        std::cout << "(" << nn[1][0] << "," << nn[1][1] << "), with phi (" << nnphi[1][0] << "," << nnphi[1][1] << "," << nnphi[1][2] << ")" << std::endl;
        std::cout << "(" << nn[2][0] << "," << nn[2][1] << "), with phi (" << nnphi[2][0] << "," << nnphi[2][1] << "," << nnphi[2][2] << ")" << std::endl;
        std::cout << "(" << nn[3][0] << "," << nn[3][1] << "), with phi (" << nnphi[3][0] << "," << nnphi[3][1] << "," << nnphi[3][2] << ")" << std::endl;
    }
    
    void ComplexLattice::SaveConfig(int step){
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
   
    void ComplexLattice::CalcGij(){
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
// Casey Berger
// Created: Apr 2, 2024
// Last edited: Apr 9, 2024
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

#include "CL.h"
#include "mathlib.h" //dot, cross
#include "cllib.h" //complex dot, complex cross (put those in here??)


namespace complexlangevin{
    //public functions
    //constructor
    CL::CL(int argc, char *argv[]){
        std::cout << "CL::Lattice" << std::endl;
        CL::GetInputs_(argc, argv);
        //CL::GenerateFilename_();//set output filename
    }
    //other public functions
    void CL::SetLength(int length){
        length_ = length;
        std::cout << "Lattice length set to " << length_ << std::endl;
        //Lattice::generateFilename_();
    }

    int CL::GetLength(){
        return length_;
    }
      /*  
    std::string CL::GetFilename(){
        return filename_;
    }
    
    CL::complex_field CL::GetField(int i, int j){
        return lattice_[i][j];
    }
 
    
    void CL::Initialize(){
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
    
    void CL::PrintLattice(){
        for (int i = 0; i < length_; i++){
            for (int j = 0; j < length_; j++){
                Lattice::printPhi_(i, j);
            }
        }
    }
    
    void CL::PrintNeighbors(int i, int j){
        std::array < Lattice::vertex, 4 > nn = Lattice::getNeighbors_(i,j);
        std::array < Lattice::field, 4 > nnphi = Lattice::getNeighborPhis_(i,j);
        std::cout << "At (i,j) = " << i << "," << j << " the neighbors are: " << std::endl;
        std::cout << "(" << nn[0][0] << "," << nn[0][1] << "), with phi (" << nnphi[0][0] << "," << nnphi[0][1] << "," << nnphi[0][2] << ")" << std::endl;
        std::cout << "(" << nn[1][0] << "," << nn[1][1] << "), with phi (" << nnphi[1][0] << "," << nnphi[1][1] << "," << nnphi[1][2] << ")" << std::endl;
        std::cout << "(" << nn[2][0] << "," << nn[2][1] << "), with phi (" << nnphi[2][0] << "," << nnphi[2][1] << "," << nnphi[2][2] << ")" << std::endl;
        std::cout << "(" << nn[3][0] << "," << nn[3][1] << "), with phi (" << nnphi[3][0] << "," << nnphi[3][1] << "," << nnphi[3][2] << ")" << std::endl;
    }
    
    void CL::SaveConfig(int step){
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
    */
    
    //private functions
    
    //reading in parameters, changing parameters, etc
    
    void CL::GetInputs_(int argc, char *argv[]){
        //read in parameters
        //consider making one testing flag that you turn on with #ifdef TESTING_MODE so you can just do if (testing){THING TO DO;}
        std::cout << "Function: CL::GetInputs_" << std::endl;
        #ifdef TESTING_MODE 
            std::cout << "Function: CL::GetInputs_" << std::endl;
        #endif
            std::string str, filename;
            const int n_params = 3; //L, nL, freq
            std::string inputs [n_params] = {"length","nL", "freq"};//read in keywords for parameters
            if (argc != 2){ //exits if input file is not given
                std::cerr << "Usage: ./executable_name input.txt"<< std::endl << "Exiting program" << std::endl;
                std::exit(10);
            }
            else{
                std::ifstream input_file(argv[1]);
                if (!input_file.is_open()){
                    std::cerr << "input file cannot be opened";
                    std::exit(10);
                }
                else{
                    int count = 0;
        #ifdef TESTING_MODE
                    std::cout << "Starting param search in file: ";
                    for (int n=0; n<n_params; n++){
                        std::cout << inputs[n] << ',';
                    }
                    std::cout << std::endl;
        #endif  
                    while (count < n_params) {
                        while (getline(input_file, str)) {
                            //search for params in input
                            size_t found = str.find(inputs[count]);
                            size_t start;
                            if (found != std::string::npos) {
                                start = str.find_last_of(' ');
                                inputs[count] = str.substr(start + 1);
                                count++;
                            }
                            else{
                                //if your inputs file doesn't have that parameter listed 
                                std::cerr << "parameter "<< inputs[count] << " not in input file.";
                                std::exit(10);
                            }
                        }
                    }
                    int length = stod(inputs[0]);
                    CL::SetLength(length);
                    int nL = stod(inputs[1]);
                    int freq = stod(inputs[2]);
        #ifdef TESTING_MODE
                    std::cout << "parameters acquired: ";
                    for (int n=0; n<n_params; n++){
                        std::cout << inputs[n] << ',';
                    }
                    std::cout << std::endl;
        #endif  
                }
            }
        }
    
    /*
    
    CL::field CL::makePhi_(){
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

    int CL::plusOne_(int i){
        //tested 5/30/2023
        if(i+1 == length_){ return 0;}
        else{return i+1;}
    }
    
    int CL::minusOne_(int i){
        //tested 5/30/2023
        if(i == 0){ return length_-1;}
        else{return i-1;}
    }
    
    
    std::array < CL::vertex, 4> CL::getNeighbors_(int i, int j){
        //tested 6/1/2023
        std::array < CL::vertex, 4> nn;
        nn[0][0] = CL::plusOne_(i);
        nn[0][1] = j;
        nn[1][0] = i;
        nn[1][1] = CL::plusOne_(j);
        nn[2][0] = CL::minusOne_(i);
        nn[2][1] = j;
        nn[3][0] = i;
        nn[3][1] = CL::minusOne_(j);
        return nn;
    }
    
    std::array < CL::field, 4 > CL::getNeighborPhis_(int i, int j){
        //tested 6/1/2023
        std::array < CL::field, 4> nnPhis;
        nnPhis[0] = CL::getPhi(CL::plusOne_(i), j);
        nnPhis[1] = CL::getPhi(i, CL::plusOne_(j));
        nnPhis[2] = CL::getPhi(CL::minusOne_(i), j);
        nnPhis[3] = CL::getPhi(i, CL::minusOne_(j));
        return nnPhis;
    }
    
    void CL::printField_(int i, int j){
        //tested 5/30/2023
        CL::field phi = CL::getPhi(i,j);
        std::cout << "At point (" << i << "," << j <<"), phi = (";
        std::cout << phi[0] << "," << phi[1] << ","<< phi[2] << ")" << std::endl;
    }

    void CL::generateFilename_(){
        std::string l_str   = std::to_string(length_);
        std::string b_str   = std::to_string(beta_);
        std::string th_str  = std::to_string(itheta_);
        std::string nt_str  = std::to_string(nTherm_);
        std::string nmc_str = std::to_string(nMC_);
        std::string frq_str = std::to_string(freq_);
        std::string fname = "nonlinearsigma_data_L_" + l_str + "_beta_" + b_str + "_itheta_" + th_str + "_ntherm_" + nt_str + "_nMC_" + nmc_str + "_freq_" + frq_str + ".csv";
        filename_ = fname;
    }
    */
    
}//end class definition
// Casey Berger
// Created: Apr 2, 2024
// Last edited: Apr 10, 2024
#include <array>
#include <vector>
#include <omp.h>
#pragma once
/*
What sort of conventions should you use here?

How to distinguish function names from variable names? 
    - variables are all lower case and have underscores between words
    - functions are all capitalized and have no spaces between words

Internal functions and variables end with "_" (or start with??)
*/

namespace complexlangevin{
    class CL {
        public:
        
        using complex_number = std::array<double, 2>;
        using complex_field = std::array<complex_number,3>;//is there a way to make this n dimensional??
        
        //constructor
        CL(int argc, char *argv[]); //is there a way to make this more flexible?
        
        //other public functions
        //set or update lattice parameters
        void SetLength(int length);
        
        //retrieve lattice parameters
        int GetLength();
        //std::string GetFilename();
        //field GetField(int i, int j);
    
        //initialize the lattice
        //void Initialize(); 
        
        //print things to the screen
        //void PrintLattice();
        
        //write things to file
        //void SaveConfig(int step);

    
        //private members -- only accessible within the class functions
        private:
        //variables
        int length_;
        //std::vector < std::vector < complex_field > > lattice_;
        int nL_;
        int freq_;
        //std::string filename_;
        
        //functions
        void GetInputs_(int argc, char *argv[]);
        //complex_field MakePhi_(); 
        //int PlusOne_(int i); 
        //int MinusOne_(int i);
        //std::array < vertex, 4 > GetNeighbors_(int i, int j);
        //std::array < complex_field, 4 > GetNeighborFields_(int i, int j);
        //void PrintField_(int i, int j); 
        //void GenerateFilename_();
    };  
}
// Casey Berger
// Created: May 24, 2023
// Last edited: May 24, 2023
// This needs to: define a phi at each point in the lattice, pick out phi values given a location on the lattice, calculate dot and cross products.
#pragma once

#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>


#ifndef LATTICE_H
#define LATTICE_H

class Lattice {
    public:
        Lattice(int length);
        void initialize();
    
    private:
        double length_;

};
#endif
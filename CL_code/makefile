SHELL=/bin/sh


# set flags for compilation

TESTING_MODE ?= FALSE
USE_OMP ?= FALSE
USE_GPROF ?= FALSE
USE_ARCCOS ?= FALSE
USE_CONST_RN ?= FALSE


# create flags for compilation

ifeq ($(TESTING_MODE), TRUE)
	DEFINES += -DTESTING_MODE
endif

ifeq ($(USE_GPROF),TRUE)
	DEFINES += -pg
endif

ifeq ($(USE_ARCCOS),TRUE)
	DEFINES += -DUSE_ARCCOS
endif

ifeq ($(USE_CONST_RN),TRUE)
	DEFINES += -DTEST_CONSTANT_RN
endif

CXX = g++

ifeq ($(USE_OMP), TRUE)
	CXXFLAGS = -Wall -O3 -g -std=c++0x -fopenmp
else
	CXXFLAGS = -Wall -g -std=c++0x --coverage
endif

CXXFLAGS += $(DEFINES)

# define command to make all

all: CL test

CL: CL_simulation.o CL.o cllib.o mathlib.o
	$(CXX) $(CXXFLAGS) -o CL CL_simulation.o CL.o cllib.o mathlib.o

CL_simulation.o: CL_simulation.cpp CL.h cllib.h mathlib.h 
	$(CXX) $(CXXFLAGS) -c CL_simulation.cpp 

test: test.o CL.o cllib.o mathlib.o
	$(CXX) $(CXXFLAGS) -o test test.o CL.o cllib.o mathlib.o

test.o: test.cpp CL.h cllib.h mathlib.h 
	$(CXX) $(CXXFLAGS) -c test.cpp 

# define dependencies
CL.o: CL.cpp CL.h cllib.h mathlib.h
	$(CXX) $(CXXFLAGS) -c CL.cpp
    
cllib.o: cllib.cpp cllib.h mathlib.h
	$(CXX) $(CXXFLAGS) -c cllib.cpp
    
mathlib.o: mathlib.cpp mathlib.h
	$(CXX) $(CXXFLAGS) -c mathlib.cpp

# define command to make clean

clean:
	rm -f *.o 
	rm -f *.gcno
	rm -f *.gcda
	touch *


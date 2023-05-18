CXX=g++


default: lattice_class.o qmc_runner.o
	$(CXX) -o run/qmc_runner lattice_class.o qmc_runner.o

lattice_class.o: lattice_class.cpp
	$(CXX) -c lattice_class.cpp

qmc_runner.o: qmc_runner.cpp
	$(CXX) -c qmc_runner.cpp

nonlinearsigma: nonlinearsigma.o lattice.o action_suite.o monte_carlo.o test_suite.o
	$(CXX) -o nonlinearsigma nonlinearsigma.o lattice.o action_suite.o monte_carlo.o test_suite.o

nonlinearsigma.o: nonlinearsigma.cpp lattice.h action_suite.h monte_carlo.h test_suite.h
	$(CXX) -c nonlinearsigma.cpp 

# Dependencies
lattice.o: lattice.cpp lattice.h 
	$(CXX) -c lattice.cpp
    
action_suite.o: action_suite.cpp action_suite.h
	$(CXX) -c action_suite.cpp
   
monte_carlo.o: monte_carlo.cpp monte_carlo.h
	$(CXX) -c monte_carlo.cpp
    
test_suite.o: test_suite.cpp test_suite.h
	$(CXX) -c test_suite.cpp

clean:
	rm *.o

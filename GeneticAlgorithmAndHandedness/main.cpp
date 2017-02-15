//
//  main.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//  Adapted from code provided by Arash


#include "GenAlgoAndHands.h"
using namespace std;

int main(int argc, const char * argv[]) {
    
    if (! (argc == 4 || argc == 5 || argc == 6)){
        std::cout << "Usage: " << argv[0] << " option (0, 1, 2, 3, or 4) density alpha [n_iterations=2000  popCount=5000]\n";
        exit(-2);
    }
   

    // write log of command that calls functino 
    ofstream runlog("run.log",std::ofstream::out);
    for (unsigned int i=0;i<argc;i++) runlog << argv[i] << " ";
    runlog << "\n";
    runlog.close();
    
    const int s_count = atoi(argv[1]);
    double density = atof(argv[2]);
    double alpha = atof(argv[3]);

    int n_iterations=2000;
    int popCount=5000;

    if (argc >= 5) n_iterations = atoi(argv[4]);
    if (argc == 6) popCount = atoi(argv[5]);
    
    assert(density > 0.);
    assert(abs(alpha)<1.0);
    assert(n_iterations > 1);
    assert(popCount > 2);
    
    switch (s_count) {
        case 0:
            TestNmolecule1(density,alpha,n_iterations,popCount);
            break;
        case 1:
            TestNmolecule2(density,alpha,n_iterations,popCount);
            break;
        case 2:
            TestNmolecule2R(density,alpha,n_iterations,popCount);
            break;
        case 3:
            TestNmolecule4R(density, alpha,n_iterations,popCount);
            break;
        case 4:
            TestNmolecule2E(density, alpha,n_iterations,popCount);
            break;
        default:
            std::cerr << "Unknown choice for switch option\n";
            break;
    }
    
    return 0;
}


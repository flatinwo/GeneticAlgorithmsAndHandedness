//
//  search_drivers.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/6/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include "search_drivers.hpp"
#include "search_driver.hpp"
#include "crystal_search.hpp"


void TestNmolecule1(double density=0.24, double alpha=0.5,
                    int iterations=2000, int popCount=5000){
    const short n=1;
    SearchDriver onemol (n,
                        SearchDriver::SEARCHTYPE::ENANTIOPURE,
                        SearchDriver::BITTYPE::RANDOMSWEEP,
                        density,
                        alpha,
                        popCount,
                        iterations,
                        7,
                        31);
    onemol.run();
    
}


void TestNmolecule2(double density=0.24, double alpha=0.5,
                    int iterations=2000, int popCount=5000){
    
    const short n=2;
    SearchDriver twomol (n,
                         SearchDriver::SEARCHTYPE::RACEMIZING,
                         SearchDriver::BITTYPE::RANDOMSWEEP,
                         density,
                         alpha,
                         popCount,
                         iterations,
                         7,
                         31);
    
    twomol.run();
}

void TestNmolecule2E(double density=0.24, double alpha=0.5,
                     int iterations=2000, int popCount=5000){
    const short n=2;
    SearchDriver twomol (n,
                         SearchDriver::SEARCHTYPE::ENANTIOPURE,
                         SearchDriver::BITTYPE::RANDOMSWEEP,
                         density,
                         alpha,
                         popCount,
                         iterations,
                         7,
                         31);
    
    twomol.run();
}


void TestNmolecule2R(double density=0.24, double alpha=0.5,
                     int iterations=2000, int popCount=5000){
    
    const short n=2;
    SearchDriver twomol (n,
                         SearchDriver::SEARCHTYPE::RACEMIC,
                         SearchDriver::BITTYPE::RANDOMSWEEP,
                         density,
                         alpha,
                         popCount,
                         iterations,
                         7,
                         31);
    
    twomol.run();
    
    
}

void TestNmolecule4R(double density=0.24, double alpha=0.5,
                     int iterations=2000, int popCount=5000){
    return;
    
    
}

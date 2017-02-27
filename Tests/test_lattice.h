//
//  test_lattice.h
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/24/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#ifndef test_lattice_h
#define test_lattice_h

#include "gene.hpp"
#include "gtest/gtest.h"
#include "GenAlgoAndHands.h"
#include <list>
#include "lattice.hpp"

// I want to write a test to see if x, y, and angles follow
// rule, basically randomize lattice vector and then check to see if tests
//.

typedef lattice Lattice;

struct LatticeTest: public testing::Test{
    Lattice* mylattice;
    std::list<double>* densities;
    LatticeTest(){
        mylattice = new Lattice();
        double mydoubles[] = {0.10,0.24,0.82,1.62};
        densities = new list<double>(mydoubles, mydoubles + sizeof(mydoubles) / sizeof(double) );
        
    }
    
    virtual ~LatticeTest(){
        delete mylattice;
    }
    
};

TEST_F(LatticeTest, densityOrhoTest){
    
    
    
    mylattice->updateVectors(0.8, 0.7, PI/2, PI, PI/2.);
    for (std::list<double>::iterator it=densities->begin(); it != densities->end(); it++){
        mylattice->setDensity(*it);
        double expected_x = (mylattice->getLatticeVector(1).getLength())/
        (mylattice->getLatticeVector(0).getLength());
        double expected_xy =(mylattice->getLatticeVector(2).getLength())/
        (mylattice->getLatticeVector(0).getLength());
        
        EXPECT_NEAR(expected_x,mylattice->getX(),0.0001);
        EXPECT_NEAR(expected_xy,
                    (mylattice->getX())*(mylattice->getY())
                    ,0.0001);

        
    }
}

TEST_F(LatticeTest, densityMonoTest){
    
    mylattice->updateVectors(0.8, 0.7, PI/2, 0.1, PI/2.);
    for (std::list<double>::iterator it=densities->begin(); it != densities->end(); it++){
        mylattice->setDensity(*it);
        double expected_x = (mylattice->getLatticeVector(1).getLength())/
        (mylattice->getLatticeVector(0).getLength());
        double expected_xy =(mylattice->getLatticeVector(2).getLength())/
        (mylattice->getLatticeVector(0).getLength());
        
        EXPECT_NEAR(expected_x,mylattice->getX(),0.0001);
        EXPECT_NEAR(expected_xy,
                    (mylattice->getX())*(mylattice->getY())
                    ,0.0001);
        
        
    }
}

#endif /* test_lattice_h */

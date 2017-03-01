//
//  molecules.h
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/28/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#ifndef molecules_h
#define molecules_h

#include "operators.hpp"

const struct AllMolecules{
    
    std::vector<double3> LTetramer, DTetramer;
    
    AllMolecules(){
        LTetramer.reserve(4);
        LTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
        LTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
        LTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
        LTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
        
        DTetramer.reserve(4);
        DTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
        DTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
        DTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
        DTetramer[3] = double3( 0.52915050,  0.26457525, 1.32287625);
        
    }
    
} moleculestore;




#endif /* molecules_h */

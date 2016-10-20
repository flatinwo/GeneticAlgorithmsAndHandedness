//
//  individual.hpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef individual_hpp
#define individual_hpp

#include <iostream>
#include <vector>
#include "gene.hpp"

using namespace std;

class individual
{
private:
    
public:
    double fitness;
    double relativeFitness;
    vector<gene> genes;
    
    // constructors
    individual();
    individual(vector<gene> _genes);
    
    // force decode by setting decoded=false
    void forceDecode();
    
    // overload operators
    bool operator < (const individual& _ind) const;
    bool operator > (const individual& _ind) const;
};

#endif /* individual_hpp */

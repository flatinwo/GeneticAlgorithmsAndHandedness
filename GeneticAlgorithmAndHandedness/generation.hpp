//
//  generation.hpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef generation_hpp
#define generation_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include "individual.hpp"
#include "ran2.hpp"

using namespace std;

class generation
{
private:
    unsigned int generationCount;
    
    // renormalize fitness, so that maximum fitness is 1
    void renormalizeFitness();
    
    // check whether suitable parents exist
    bool checkParents();
    
    // swap fittest individual to first index
    void promoteFittest();
public:
    bool elitism; // if elitism is set to true, then the fittest individual (idx 0) is left untouched
    vector<individual> individuals;
    
    // constructors
    generation();
    generation(vector<individual> _individuals);
    
    // functions for GA
    // calculate fitness
    void calcFitness(double(*f)(individual));
    
    // randomizes genetic data;
    // _pIndividual is probability to pick an individual
    // _pGene is probability to randomize a gene
    void randomize(double _pIndividual=1.0, double _pGene=1.0);
    
    // perform point mutation;
    // _pIndividual is probability to pick an individual
    // _pFlip is probability to flip a bit in the code
    void mutate(double _pIndividual, double _pFlip);
    
    // mate individuals using the fitness proportional selection
    void mateFitnessProportional();
    
    // sort individuals by fitness
    void sortByFitness();
};

#endif /* generation_hpp */

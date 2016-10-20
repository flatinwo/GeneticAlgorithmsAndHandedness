//
//  individual.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "individual.hpp"

individual::individual()
{
    fitness = 0.0;
}

individual::individual(vector<gene> _genes)
{
    genes = _genes;
    fitness = 0.0;
}

void individual::forceDecode()
{
    for (unsigned int i=0; i<genes.size(); i++)
        genes[i].forceDecode();
}

bool individual::operator < (const individual& _ind) const
{
    return (fitness < _ind.fitness);
}

bool individual::operator > (const individual& _ind) const
{
    return (fitness > _ind.fitness);
}

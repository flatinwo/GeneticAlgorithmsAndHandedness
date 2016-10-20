//
//  generation.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "generation.hpp"

generation::generation()
{
    generationCount = 0;
    elitism = false;
}

generation::generation(vector<individual> _individuals)
{
    generationCount = 0;
    elitism = false;
    individuals = _individuals;
}


// renormalize fitness, so that maximum fitness is 1
void generation::renormalizeFitness()
{
    double fitnessSum = 0.0;
    for (unsigned int i=0; i<individuals.size(); i++)
    {
        fitnessSum += individuals[i].fitness;
    }
    
    if (fitnessSum > 0.0)
    {
        for (unsigned int i=0; i<individuals.size(); i++)
        {
            individuals[i].relativeFitness = individuals[i].fitness/fitnessSum;
        }
    }
}

bool generation::checkParents()
{
    renormalizeFitness();
    
    const double fitnessThreshold = 1E-6;
    unsigned int nonZeroFitnessCount = 0;
    
    for (unsigned int i=0; i<individuals.size(); i++)
    {
        if (individuals[i].relativeFitness > fitnessThreshold)
        {
            nonZeroFitnessCount++;
        }
    }
    
    if (nonZeroFitnessCount < 2)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void generation::promoteFittest()
{
    double fitnessMax = 0.0;
    unsigned int fitnessMaxIndex = 0;
    for (unsigned int i=0; i<individuals.size(); i++)
    {
        if (individuals[i].fitness > fitnessMax)
        {
            fitnessMaxIndex = i;
            fitnessMax = individuals[i].fitness;
        }
    }
    individual tmp = individuals[0];
    individuals[0] = individuals[fitnessMaxIndex];
    individuals[fitnessMaxIndex] = tmp;
}

void generation::randomize(double _pIndividual, double _pGene)
{
    for (unsigned int i=0; i<individuals.size(); i++)
    {
        if (ran2(&seed) <= _pIndividual)
        {
            for (unsigned int j=0; j<individuals[i].genes.size(); j++)
            {
                if (ran2(&seed) <= _pGene)
                {
                    individuals[i].genes[j].randomize();
                }
            }
        }
    }
}

// Calculate fitness of all the individuals in this generation
// We use the externally provided function f here
void generation::calcFitness(double(*f)(individual _ind))
{
    for (unsigned int i=0; i<individuals.size(); i++)
    {
        individuals[i].fitness = (*f)(individuals[i]);
    }
    renormalizeFitness();
}

// Perform point mutation of all individuals in the generation
// Skip the fittest one (at index 0) if elitism is enabled
void generation::mutate(double _pIndividual, double _pFlip)
{
    unsigned int startIndex = 0;
    
    if (elitism)
    {
        startIndex = 1;
        promoteFittest();
    }
    
    for (unsigned int i=startIndex; i<individuals.size(); i++)
    {
        if (ran2(&seed) <= _pIndividual)
        {
            for (unsigned int j=0; j<individuals[i].genes.size(); j++)
            {
                // flip bits
                individuals[i].genes[j].flip(_pFlip);
            }
        }
    }
}

void generation::mateFitnessProportional()
{
    // first check whether suitable number of mates (2) exist with fitness > 0
    if (!checkParents())
    {
        cerr<<"Warning: No suitable parents found in generation "<<generationCount<<". Increase population or make fitness function softer."<<endl;
        return;
    }
    
    // prepare vectors for children
    unsigned int childCount = individuals.size();
    
    if (elitism)
    {
        childCount--;
        promoteFittest();
    }
    
    vector<individual> children;
    children.reserve(childCount);
    
    // calculate relative fitness
    renormalizeFitness();
    
    for (unsigned int i=0; i<childCount; i++)
    {
        bool fatherFound = false;
        bool motherFound = false;
        unsigned int fatherIndex = 0;
        unsigned int motherIndex = 0;
        
        do
        {
            for (unsigned int j=0; j<individuals.size(); j++)
            {
                if (ran2(&seed) < individuals[j].relativeFitness)
                {
                    fatherFound = true;
                    fatherIndex = j;
                    break;
                }
            }
        } while (!fatherFound);
        
        do
        {
            for (unsigned int j=0; j<individuals.size(); j++)
            {
                if ((ran2(&seed) < individuals[j].relativeFitness) && (j != fatherIndex))
                {
                    motherFound = true;
                    motherIndex = j;
                    break;
                }
            }
        } while (!motherFound);
        
        // create child by combining genetic data of father and mother
        // determine intersection point
        unsigned int intersectGene = ran2(&seed)*individuals[fatherIndex].genes.size();
        unsigned int intersectCode = ran2(&seed)*individuals[fatherIndex].genes[intersectGene].code.size();
        
        // first create with father's genetic data;
        individual child = individuals[fatherIndex];
        child.fitness = 0.0;
        child.relativeFitness = 0.0;
        child.forceDecode();
        
        // now replace with mother's genetic data from intersection point onwards
        for (unsigned int j=intersectGene; j<child.genes.size(); j++)
        {
            unsigned int startIndex = 0; // at intersected gene, start at intersected code; otherwise use full code
            if (j == intersectGene)
            {
                startIndex = intersectCode;
            }
            
            for (unsigned int k=startIndex; k<child.genes[j].code.size(); k++)
            {
                child.genes[j].code[k] = individuals[motherIndex].genes[j].code[k];
            }
        }
        
        children.push_back(child);
    }
    
    // replace old generation with new generation, and take into account elitisim
    for (unsigned int i=0; i<childCount; i++)
    {
        individuals[individuals.size()-childCount+i] = children[i];
    }
    
    generationCount++;
}

void generation::sortByFitness()
{
    sort(individuals.rbegin(), individuals.rend()); // descending order ; otherwise use .begin() and .end()
}

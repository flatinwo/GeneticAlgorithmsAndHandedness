//
//  gene.hpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef gene_hpp
#define gene_hpp

#include <vector>
#include <iostream>
#include <string>
#include "ran2.hpp"

using namespace std;

class gene
{
private:
    unsigned int length; // length of the code, >1
    double maxIntValue; // maximum integer value the gene can take using binary code
    double minValue; // minimum real value the gene should take
    double maxValue; // maximum real vlaue the gene should take
    double value; // real value
    bool decoded; // flag that indicates whether decoded value has been stored or not
    
    void initGene(unsigned int _length, double _minValue, double _maxValue);
    void decode();
    void encode(double _value);
public:
    vector<bool> code; // this is the actual code, consisting of 0s and 1s
    
    // constructors
    gene(unsigned int _length, double _minValue, double _maxValue);
    gene(unsigned int _length, double _value, double _minValue, double _maxValue);
    
    // get and set functions
    double getValue();
    bool setValue(double);
    
    // force decode by setting decoded=false
    void forceDecode();
    
    // randomize genetic code
    void randomize();
    
    // flip parts of genetic code
    void flip(double _pFlip);
    
    // return code as a string
    string makeString();
};

#endif /* gene_hpp */

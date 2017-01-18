//
//  gene.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "gene.hpp"

gene::gene(unsigned int _length, double _minValue, double _maxValue)
{
    initGene(_length, _minValue, _maxValue);
    encode(_minValue);
}

gene::gene(unsigned int _length, double _value, double _minValue, double _maxValue)
{
    initGene(_length, _minValue, _maxValue);
    encode(_value);
}

void gene::initGene(unsigned int _length, double _minValue, double _maxValue)
{
    if ((_length > 0) && (_length < 32))
    {
        length = _length;
    }
    else if (_length >= 32)
    {
        length = 31;
        cerr<<"Warning: Maximum length per gene is 32"<<endl;
    }
    else
    {
        length = 1;
        cerr<<"Warning: Minimum length per gene is 1"<<endl;
    }
    
    minValue = _minValue;
    maxValue = _maxValue;
    code.assign(length, 0);
    
    // determine number of different values, the code can hold
    maxIntValue = 0;
    for (unsigned int i=0; i<length; i++)
    {
        maxIntValue += 1 << i;
    }
}

// decode binary code to a value in the range [minValue,maxValue] with stepsize 1/maxIntValue
void gene::decode()
{
    double r = 0.0;
    for (unsigned int i=0; i<length; i++) // interpret binary code as an unsigned integer
    {
        if (code[i])
        {
            r += 1 << i;
        }
    }
    r /= maxIntValue; // normalize code and transform into desired range
    value = minValue + r*(maxValue-minValue);
    
    decoded = true;
}

void gene::encode(double _value)
{
    // first, map value to range between 0...maxIntValue
    double r = maxIntValue*(_value - minValue)/(maxValue-minValue);
    double remainder = r;
    
    for (unsigned int i=0; i<length; i++)
    {
        const unsigned int idx = length-i-1;
        double divValue = 1 << idx; // determine value of a binary number with 1 at position idx
        double fraction = remainder/divValue; // this value should be < 2
        unsigned int fractionIndex = fraction; // this should be either 1 or 0
        remainder = (fraction-fractionIndex)*divValue;
        
        code[idx] = fractionIndex;
    }
    
    decoded = false;
}

double gene::getValue()
{
    if (decoded)
    {
        return value;
    }
    else
    {
        decode();
        return value;
    }
}

bool gene::setValue(double _value)
{
    if ((_value < minValue) || (_value > maxValue))
    {
        return false; // value to be set is out of range
    }
    else
    {
        encode(_value);
        return true;
    }
}

void gene::forceDecode()
{
    decoded = false;
}

void gene::randomize()
{
    for (unsigned int i=0; i<length; i++)
    {
        if (ran2(&seed) < 0.5)
        {
            code[i] = 0;
        }
        else
        {
            code[i] = 1;
        }
    }
    decoded = false;
}

void gene::flip(double _pFlip)
{
    for (unsigned int i=0; i<length; i++)
    {
        if (ran2(&seed) < _pFlip)
        {
            code[i] = !code[i];
        }
    }
    decoded = false;
}

string gene::makeString()
{
    string codeString = "";
    for (unsigned int i=0; i<length; i++)
    {
        if (code[i])
        {
            codeString += "1";
        }
        else
        {
            codeString += "0";
        }
    }
    return codeString;
}

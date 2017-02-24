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

    updateLength(_length);
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

void gene::updateLength(unsigned int _length){
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
}

void gene::updateBittage(unsigned int _length){
    double curr_val = getValue();
    initGene(_length, minValue, maxValue);
    encode(curr_val);
}

// decode binary code to a value in the range [minValue,maxValue] with stepsize 1/maxIntValue
void gene::decode(bool diagnoze)
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
    
    if (!diagnoze) gval.value = value;
    
    decoded = true;
}

void gene::encode(double _value, bool diagnoze)
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
    
    if (!diagnoze){
        gval.code = code;
        gval.codeString = makeString();
    }
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

void gene::refreshToOldState(){
    code = gval.code;
    codeString = gval.codeString;
    value = gval.value;
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

double gene::returnValAndrefresh(){
    decode(true);
    double val = getValue();
    forceDecode();
    refreshToOldState();
    return val;
}

// gives resolution of gene based on bit
double gene::getResolution(){
    code.assign(code.size(),false);
    code.front() = true;
    return returnValAndrefresh();

}// provide me with the resolution


string gene::getResolutionInBits(){
    //encode(getResolution(),true);
    code.assign(code.size(),false); //0.00879765
    code.front() = true;//calling get resolution will give rounding errors
    tempstr = makeString();
    refreshToOldState();
    return tempstr;
}

// self descriptive
double gene::getMaxValue(){
    code.assign(code.size(),true);
    return returnValAndrefresh();
}

string gene::getMaxValueInBits(){
    encode(maxValue,true);
    tempstr = makeString();
    forceDecode();
    refreshToOldState();
    return tempstr;

}

// self descriptive
double gene::getMinValue(){
    code.assign(code.size(),false);
    return returnValAndrefresh();
}

string gene::getMinValueInBits(){
    encode(minValue,true);
    tempstr = makeString();
    forceDecode();
    refreshToOldState();
    return tempstr;
}

string gene::makeString()
{
    codeString = "";
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

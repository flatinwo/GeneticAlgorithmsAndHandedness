//
//  operators.hpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 12/12/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef operators_hpp
#define operators_hpp

#include <stdio.h>
#include "lattice.hpp"
#include <vector>
#include <map>
#include <string>

using namespace std;

double4 axisangle(const double4& quat);
vector<double3> rotate(const double4& quat, const vector<double3>& points);
void zeroCOM(vector<double3>& mol);
void normalize(double& q0, double& q1, double& q2, double& q3);
void normalizequaternion(double& q0, double& q1, double& q2, double& q3);

/**
 @param ratio2 should be (sigma/distance)^2
 @return Lennard-Jones associated with ratio, assumes epsilon=1.
 */
double vLJ(const double& ratio2);

typedef std::pair< vector<double3>, std::string > molecule;

template<unsigned int n>
struct factorial {
    enum { value = n * factorial<n-1>::value };
};

template<>
struct factorial<0>{
    enum {value = 1};
};

template <unsigned int n, unsigned int k>
struct combination{
    enum { value = factorial<n>::value/
        (factorial<n-k>::value*factorial<k>::value) };
};


#endif /* operators_hpp */

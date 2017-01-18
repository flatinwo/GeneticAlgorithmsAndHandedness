//
//  operators.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 12/12/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "operators.hpp"
#include <vector>
#include <cassert>

using namespace std;

double vLJ(const double& ratio2){
    double ratio6 = pow(ratio2,3.0);
    return 4.0*ratio6*(ratio6 - 1.0);
}

double4 axisangle(const double4& quat){
    double x,y,z;
    double angle = 2.*acos(quat.scalar);
    double s = sqrt(1.0-quat.scalar*quat.scalar);
    if (s < 0.001) { // test to avoid divide by zero, s is always positive due to sqrt
        // if s close to zero then direction of axis not important
        x = quat.coords.x; // if it is important that axis is normalised then replace with x=1; y=z=0;
        y = quat.coords.y;
        z = quat.coords.z;
    } else {
        x = quat.coords.x / s; // normalise axis
        y = quat.coords.y / s;
        z = quat.coords.z / s;
    }
    
    return double4(x,y,z,angle);
}

vector<double3> rotate(const double4& quat, const vector<double3>& points){
    vector<double3> result(points.size(),double3(0.,0.,0.));
    double a,b,c,d;
    a = quat.scalar;
    b = quat.coords.x;
    c = quat.coords.y;
    d = quat.coords.z;
    
    double rotationMatrix[3][3]={0.};
    
    
    rotationMatrix[0][0] = a*a + b*b - c*c - d*d; rotationMatrix[0][1] = 2*(b*c - a*d);
    rotationMatrix[0][2] = 2.*(b*d + a*c);
    
    
    rotationMatrix[1][0] = 2.*(b*c + a*d); rotationMatrix[1][1] = a*a - b*b + c*c - d*d;
    rotationMatrix[1][2] = 2.*(c*d - a*b);
    
    rotationMatrix[2][0] = 2.*(b*d - a*c); rotationMatrix[2][1] = 2.*(c*d + a*b);
    rotationMatrix[2][2] = a*a - b*b - c*c + d*d;
    
    for (unsigned int i=0; i<result.size(); i++) {
        result[i].x = rotationMatrix[0][0]*points[i].x + rotationMatrix[0][1]*points[i].y + rotationMatrix[0][2]*points[i].z;
        result[i].y = rotationMatrix[1][0]*points[i].x + rotationMatrix[1][1]*points[i].y + rotationMatrix[1][2]*points[i].z;
        result[i].z = rotationMatrix[2][0]*points[i].x + rotationMatrix[2][1]*points[i].y + rotationMatrix[2][2]*points[i].z;
    }
    return result;
}

void zeroCOM(vector<double3>& mol){
    double3 com(0.,0.,0.);
    for (auto& i: mol) com = com + i;
    com.x /= 4; com.y /= 4; com.z /= 4;
    for (auto& i : mol) i = i - com;
}


void normalize(double& q0, double& q1, double& q2, double& q3){
    double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    if (norm > 0.001) {
        q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
    }
    else{
        q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
        assert( sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3) > abs(1.-0.0001) );
        //std::cerr << "Attempting to normalize zero vector\n";
        //exit(-1);
    }

}

void normalizequaternion(double& q0, double& q1, double& q2, double& q3){
    q3 = sqrt(1. - (q0*q0 + q1*q1 + q2*q2));
}


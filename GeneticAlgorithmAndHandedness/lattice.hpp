//
//  lattice.hpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef PI
#define PI 3.141592653589793238462643
#endif

#ifndef lattice_hpp
#define lattice_hpp

#include <iostream>
#include <cmath>

using namespace std;

struct double3
{
    double x, y, z;
    double3() : x(0.0), y(0.0), z(0.0) {}
    double3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    
    double3 operator+(const double3& b)	{return double3(this->x + b.x, this->y + b.y, this->z + b.z);} // Sum
    double3 operator-(const double3& b)	{return double3(this->x - b.x, this->y - b.y, this->z - b.z);} // Difference
    double operator*(const double3& b) {return (this->x*b.x + this->y*b.y + this->z*b.z);} // Scalar product
    double3 operator%(const double3& b) // Cross product
    {
        return double3(this->y*b.z - this->z*b.y, this->z*b.x - this->x*b.z, this->x*b.y - this->y*b.x);
    }
    double get2Norm(){ return (x*x + y*y + z*z); }
    double getLength(){ return sqrt(get2Norm());}
};

struct double4{
    double3 coords;
    double scalar;
    
    double4():coords(0,0,0),scalar(0.){};
    double4(double _x,double _y,double _z,double _s):coords(_x,_y,_z),scalar(_s){};
    double4(double3 _coords, double _scalar=0.):coords(_coords),scalar(_scalar){};
    
};

struct lattice
{
private:
    double x, y, phi, psi, theta; // These parameters fully describe the lattice
    double a; // The lattice constant just depends on the density
    double3 a1, a2, a3; // These are the three lattice vectors transformed in Cartesian coordinates
    
    void minimizeArea();
    void sortLatticeVectors();
public:
    // Constructors
    lattice();
    lattice(double, double ,double, double, double);
    
    // Functions for processing and updating lattice
    void updateVectors(double, double, double, double, double); // Update vectors based on x, y, phi, psi, theta
    double getVolume();
    double getArea();
    double getArea(double3, double3, double3);
    double3 getLatticeVector(unsigned int);
    void setDensity(double);
    
    double getX(){ return x;};
    double getY(){ return y;};
    double getPhi(){return phi;};
    double getPsi(){return psi;};
    double getTheta(){return theta;};
};

#endif /* lattice_hpp */

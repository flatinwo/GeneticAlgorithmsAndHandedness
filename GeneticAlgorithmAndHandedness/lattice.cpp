//
//  lattice.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "lattice.hpp"

lattice::lattice()
{
    a1 = double3(1.0, 0.0, 0.0);
    a2 = double3(0.0, 1.0, 0.0);
    a3 = double3(0.0, 0.0, 1.0);
    setDensity(1.0);
}

lattice::lattice(double _x, double _y, double _phi, double _psi, double _theta)
{
    updateVectors(_x, _y, _phi, _psi, _theta);
    setDensity(1.0);
}

void lattice::updateVectors(double _x, double _y, double _phi, double _psi, double _theta)
{
    x = _x;
    y = _y;
    phi = _phi;
    psi = _psi;
    theta = _theta;
    
    if ((x != 0.0) && (y != 0.0) && (phi != 0.0) && (psi != 0.0) && (theta != 0.0))
    {
        a1 = double3(1.0, 0.0, 0.0);
        a2 = double3(x*cos(phi), x*sin(phi), 0.0);
        a3 = double3(x*y*cos(psi)*cos(theta), x*y*sin(psi)*cos(theta), x*y*sin(theta));
    }
    else
    {
        a1 = double3(1.0, 0.0, 0.0);
        a2 = double3(0.0, 1.0, 0.0);
        a3 = double3(0.0, 0.0, 1.0);
    }
    
    minimizeArea();
    sortLatticeVectors();
}

double lattice::getArea() // Calculate the are of the 3 lattice vectors
{
    return getArea(a1, a2, a3);
}

double lattice::getArea(double3 _a1, double3 _a2, double3 _a3)
{
    const double3 p1 = _a1%_a2;
    const double3 p2 = _a1%_a3;
    const double3 p3 = _a2%_a3;
    const double p1l = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
    const double p2l = sqrt(p2.x*p2.x + p2.y*p2.y + p2.z*p2.z);
    const double p3l = sqrt(p3.x*p3.x + p3.y*p3.y + p3.z*p3.z);
    
    return (p1l + p2l + p3l);
}

double lattice::getVolume()
{
    return fabs((a1%a2)*a3);
}

double3 lattice::getLatticeVector(unsigned int _index)
{
    switch (_index)
    {
        case 0:
            return double3(a*a1.x, a*a1.y, a*a1.z);
            break;
        case 1:
            return double3(a*a2.x, a*a2.y, a*a2.z);
            break;
        case 2:
            return double3(a*a3.x, a*a3.y, a*a3.z);
            break;
        default:
            return double3(a*a1.x, a*a1.y, a*a1.z);
            break;
    }
}

void lattice::setDensity(double _density)
{
    const double V = getVolume();
    a = pow(_density*V, -1.0/3.0);
}

void lattice::minimizeArea()
{
    // Calculate Area of these 12 permutations:
    // {a1+-a2, a2, a3}, {a1+-a3, a2, a3}
    // {a1, a2+-a1, a3}, {a1, a2+-a3, a3}
    // {a1, a2, a3+-a1}, {a1, a2, a3+-a2}
    // And repeat this procedure until area raises
    
    bool repeat = true;
    double m = 1024.0;
   	
    // The minimization is carried out using a binary search, where the stepping is halved after each failed minimization step
    do
    {
        const double oldArea = getArea();
       	double bestArea = oldArea;
        double newArea;
        double3 swapVector; // This holds the permutated trial vector
        double3 bestVector; // This holds the best permutated trial vector
        unsigned int swapVectorIndex = 0; // This variable indicates, whether a1, a2 or a3 has been permutated
        
        // {a1+a2, a2, a3}
        swapVector = double3(a1.x+m*a2.x, a1.y+m*a2.y, a1.z+m*a2.z);
        newArea = getArea(swapVector, a2, a3);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 0;
            bestVector = swapVector;
        }
        
        // {a1-a2, a2, a3}
        swapVector = double3(a1.x-m*a2.x, a1.y-m*a2.y, a1.z-m*a2.z);
        newArea = getArea(swapVector, a2, a3);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 0;
            bestVector = swapVector;
        }
        
        // {a1+a3, a2, a3}
        swapVector = double3(a1.x+m*a3.x, a1.y+m*a3.y, a1.z+m*a3.z);
        newArea = getArea(swapVector, a2, a3);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 0;
            bestVector = swapVector;
        }
        
        // {a1-a3, a2, a3}
        swapVector = double3(a1.x-m*a3.x, a1.y-m*a3.y, a1.z-m*a3.z);
        newArea = getArea(swapVector, a2, a3);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 0;
            bestVector = swapVector;
        }
        
        // {a1, a2+a1, a3}
        swapVector = double3(a2.x+m*a1.x, a2.y+m*a1.y, a2.z+m*a1.z);
        newArea = getArea(a1, swapVector, a3);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 1;
            bestVector = swapVector;
        }
        
        // {a1, a2-a1, a3}
        swapVector = double3(a2.x-m*a1.x, a2.y-m*a1.y, a2.z-m*a1.z);
        newArea = getArea(a1, swapVector, a3);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 1;
            bestVector = swapVector;
        }
        
        // {a1, a2+a3, a3}
        swapVector = double3(a2.x+m*a3.x, a2.y+m*a3.y, a2.z+m*a3.z);
        newArea = getArea(a1, swapVector, a3);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 1;
            bestVector = swapVector;
        }
        
        // {a1, a2-a3, a3}
        swapVector = double3(a2.x-m*a3.x, a2.y-m*a3.y, a2.z-m*a3.z);
        newArea = getArea(a1, swapVector, a3);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 1;
            bestVector = swapVector;
        }
        
        // {a1, a2, a3+a1}
        swapVector = double3(a3.x+m*a1.x, a3.y+m*a1.y, a3.z+m*a1.z);
        newArea = getArea(a1, a2, swapVector);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 2;
            bestVector = swapVector;
        }
        
        // {a1, a2, a3-a1}
        swapVector = double3(a3.x-m*a1.x, a3.y-m*a1.y, a3.z-m*a1.z);
        newArea = getArea(a1, a2, swapVector);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 2;
            bestVector = swapVector;
        }
        
        // {a1, a2, a3+a2}
        swapVector = double3(a3.x+m*a2.x, a3.y+m*a2.y, a3.z+m*a2.z);
        newArea = getArea(a1, a2, swapVector);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 2;
            bestVector = swapVector;
        }
        
        // {a1, a2, a3-a2}
        swapVector = double3(a3.x-m*a2.x, a3.y-m*a2.y, a3.z-m*a2.z);
        newArea = getArea(a1, a2, swapVector);
        if (newArea < bestArea)
        {
            bestArea = newArea;
            swapVectorIndex = 2;
            bestVector = swapVector;
        }
        
        // Use best permutation
        if (bestArea < oldArea)
        {
            switch(swapVectorIndex)
            {
                case 0:
                    a1 = bestVector;
                    break;
                case 1:
                    a2 = bestVector;
                    break;
                case 2:
                    a3 = bestVector;
                    break;
            }
        }
        else
        {
            if (m <= 1.0)
            {
                repeat = false;
            }
            else
            {
                m *= 0.5;
            }
            
            /*cout<<"Minimize area, m: "<<m<<", oldArea: "<<oldArea<<", bestArea: "<<bestArea<<endl;
             cout<<a1.x<<", "<<a1.y<<", "<<a1.z<<endl;
             cout<<a2.x<<", "<<a2.y<<", "<<a2.z<<endl;
             cout<<a3.x<<", "<<a3.y<<", "<<a3.z<<endl;*/
        }
        
    } while(repeat);
}

// Triangulate matrix a1a2a3
void lattice::sortLatticeVectors()
{
    double matrix[3][3] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    // Sort by length & insert transposed into matrix
    const double a1l = a1.x*a1.x + a1.y*a1.y + a1.z*a1.z;
    const double a2l = a2.x*a2.x + a2.y*a2.y + a2.z*a2.z;
    const double a3l = a3.x*a3.x + a3.y*a3.y + a3.z*a3.z;
    
    if ((a1l >= a2l) && (a1l >= a3l))
    {
        matrix[0][0] = a1.x; matrix[1][0] = a1.y; matrix[2][0] = a1.z;
        
        if (a2l > a3l)
        {
            matrix[0][1] = a2.x; matrix[1][1] = a2.y; matrix[2][1] = a2.z;
            matrix[0][2] = a3.x; matrix[1][2] = a3.y; matrix[2][2] = a3.z;
        }
        else
        {
            matrix[0][1] = a3.x; matrix[1][1] = a3.y; matrix[2][1] = a3.z;
            matrix[0][2] = a2.x; matrix[1][2] = a2.y; matrix[2][2] = a2.z;
        }
    }
    else if ((a2l >= a1l) && (a2l >= a3l))
    {
        matrix[0][0] = a2.x; matrix[1][0] = a2.y; matrix[2][0] = a2.z;
        
        if (a1l > a3l)
        {
            matrix[0][1] = a1.x; matrix[1][1] = a1.y; matrix[2][1] = a1.z;
            matrix[0][2] = a3.x; matrix[1][2] = a3.y; matrix[2][2] = a3.z;
        }
        else
        {
            matrix[0][1] = a3.x; matrix[1][1] = a3.y; matrix[2][1] = a3.z;
            matrix[0][2] = a1.x; matrix[1][2] = a1.y; matrix[2][2] = a1.z;
        }
    }
    else if ((a3l >= a1l) && (a3l >= a2l))
    {
        matrix[0][0] = a3.x; matrix[1][0] = a3.y; matrix[2][0] = a3.z;
        
        if (a1l > a2l)
        {
            matrix[0][1] = a1.x; matrix[1][1] = a1.y; matrix[2][1] = a1.z;
            matrix[0][2] = a2.x; matrix[1][2] = a2.y; matrix[2][2] = a2.z;
        }
        else
        {
            matrix[0][1] = a2.x; matrix[1][1] = a2.y; matrix[2][1] = a2.z;
            matrix[0][2] = a1.x; matrix[1][2] = a1.y; matrix[2][2] = a1.z;
        }
    }
    
    // Now perform QR decomposition to diagonalize matrix
    double d[3];
    
    for(unsigned int j=0; j<3; j++)
    {
        double sigma = 0.0;
        for (unsigned int i=j; i<3; i++)
        {
            sigma += matrix[i][j]*matrix[i][j];
        }
        
        double s = 0.0;
        if (matrix[j][j] < 0)
        {
            d[j] = s = sqrt(sigma);
        }
        else
        {
            d[j] = s = -sqrt(sigma);
        }
        
        double beta = 1.0/(s*matrix[j][j] - sigma);
        matrix[j][j] -= s;
        
        for(unsigned int k=j+1; k<3; k++)
        {
            double sum = 0.0;
            for (unsigned int i=j; i<3; i++)
            {
                sum += matrix[i][j]*matrix[i][k];
            }
            
            sum *= beta;
            
            for (unsigned int i=j; i<3; i++)
            {
                matrix[i][k] += matrix[i][j]*sum;
            }
        }
    }
    
    // Update lattice vectors
    a1 = double3(d[0], 0.0, 0.0);
    a2 = double3(matrix[0][1], d[1], 0.0);
    a3 = double3(matrix[0][2], matrix[1][2], d[2]);
    
    /*cout<<"Diagonalize, area: "<<getArea()<<endl;
     cout<<a1.x<<", "<<a1.y<<", "<<a1.z<<endl;
     cout<<a2.x<<", "<<a2.y<<", "<<a2.z<<endl;
     cout<<a3.x<<", "<<a3.y<<", "<<a3.z<<endl;*/
}

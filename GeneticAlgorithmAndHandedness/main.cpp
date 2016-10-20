//
//  main.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 8/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//  Adapted from code provided by Arash

#include <iostream>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "lattice.hpp"
#include "gene.hpp"
#include "individual.hpp"
#include "generation.hpp"

using namespace std;

// Define system parameters
double density = 0.24;
const double epsilon = 1.0;
const double sigma = 1.0;
const double sigma2 = sigma*sigma;
const double rCut = 6.0;
const double rCut2 = rCut*rCut;
vector< vector<double> > rotationMatrix(3,vector<double>(3,0.));
vector<double3> refTetramer(4,double3());
vector<double3> refTetramer2(refTetramer);
double3 com = double3();
double alpha = -0.5;

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

void writexyz(individual _ind){
    // Create a lattice from the supplied individual
    lattice ltc = lattice(_ind.genes[0].getValue(), _ind.genes[1].getValue(), _ind.genes[2].getValue(), _ind.genes[3].getValue(), _ind.genes[4].getValue());
    ltc.setDensity(density/2.);
    
    double dx = _ind.genes[9].getValue();;
    double dy = _ind.genes[10].getValue();
    double dz = _ind.genes[11].getValue();

    
    const double3 a1 = ltc.getLatticeVector(0);
    const double3 a2 = ltc.getLatticeVector(1);
    const double3 a3 = ltc.getLatticeVector(2);
    const double a1l = sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
    const double a2l = sqrt(a2.x*a2.x + a2.y*a2.y + a2.z*a2.z);
    const double a3l = sqrt(a3.x*a3.x + a3.y*a3.y + a3.z*a3.z);
    const double3 d1 = double3(a1.x*dx, a2.y*dy, a3.z*dz);

    
    double q0,q1,q2,q3;
    q0=_ind.genes[5].getValue();
    q1=_ind.genes[6].getValue();
    q2=_ind.genes[7].getValue();
    q3=_ind.genes[8].getValue();
    
    double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
    
    vector<double3> newTetramer = rotate(double4(q0,q1,q2,q3),refTetramer);
    
    double q00, q01, q02,q03;
    q00=_ind.genes[12].getValue();
    q01=_ind.genes[13].getValue();
    q02=_ind.genes[14].getValue();
    q03=_ind.genes[15].getValue();
    
    norm = sqrt(q00*q00 + q01*q01 + q02*q02 + q03*q03);
    if (norm > 0.001) q00 /= norm; q01 /= norm; q02 /= norm; q03 /= norm;
    
    vector<double3> newTetramer2 = rotate(double4(q00,q01,q02,q03),refTetramer2);
    for (auto& t : newTetramer2) t = t + d1;
    
    const int iMax = ceil(rCut/a1l) + 1;
    const int jMax = ceil(rCut/a2l) + 1;
    const int kMax = ceil(rCut/a3l) + 1;
    
    ofstream myxyz("lattice.xyz");
    int tetramercounted = 0;
    int ncount=0;
    
    cout << "Box info\n";
    cout << 2.*rCut/a1l << "\t" << 2.*rCut/a2l << "\t" << 2.*rCut/a3l << "\n";
    cout << 2*iMax + 1 << "\t" << 2*jMax + 1 << "\t" << 2*kMax + 1 << "\n";
    
    for (int i=-iMax; i<=iMax; i++)
    {
        for (int j=-jMax; j<=jMax; j++)
        {
            for (int k=-kMax; k<=kMax; k++)
            {
                const double3 site_ijk = double3(i*a1.x + j*a2.x + k*a3.x, i*a1.y + j*a2.y + k*a3.y, i*a1.z + j*a2.z + k*a3.z);
                if ((i==j) &&(j==k) && (k==0)){
                    cout << "I am up to " << ncount << "\n";
                }
                for (unsigned int l=0; l<newTetramer.size(); l++) {
                    double3 ijk = newTetramer[l] + site_ijk;
                    myxyz << "O\t" << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                    ncount +=3;
                }
                
                for (unsigned int l=0; l<newTetramer2.size(); l++) {
                    double3 ijk = newTetramer2[l] + site_ijk; // + d1;
                    myxyz << "N\t" << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                    ncount +=3;
                }
                tetramercounted+=2;

            }
        }
    }
    myxyz.close();

}

double fitness(individual _ind)
{
    // Create a lattice from the supplied individual
    lattice ltc = lattice(_ind.genes[0].getValue(), _ind.genes[1].getValue(), _ind.genes[2].getValue(), _ind.genes[3].getValue(), _ind.genes[4].getValue());
    
    ltc.setDensity(density/2.);
    
    double dx = _ind.genes[9].getValue();;
    double dy = _ind.genes[10].getValue();
    double dz = _ind.genes[11].getValue();

    
    const double3 a1 = ltc.getLatticeVector(0);
    const double3 a2 = ltc.getLatticeVector(1);
    const double3 a3 = ltc.getLatticeVector(2);
    const double a1l = sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
    const double a2l = sqrt(a2.x*a2.x + a2.y*a2.y + a2.z*a2.z);
    const double a3l = sqrt(a3.x*a3.x + a3.y*a3.y + a3.z*a3.z);
    const double3 d1 = double3(a1.x*dx, a2.y*dy, a3.z*dz);
    
    double q0,q1,q2,q3;
    q0=_ind.genes[5].getValue();
    q1=_ind.genes[6].getValue();
    q2=_ind.genes[7].getValue();
    q3=_ind.genes[8].getValue();
    
    double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
    
    vector<double3> newTetramer = rotate(double4(q0,q1,q2,q3),refTetramer);
    
    double q00, q01, q02,q03;
    q00=_ind.genes[12].getValue();
    q01=_ind.genes[13].getValue();
    q02=_ind.genes[14].getValue();
    q03=_ind.genes[15].getValue();
    
    norm = sqrt(q00*q00 + q01*q01 + q02*q02 + q03*q03);
    if (norm > 0.001) q00 /= norm; q01 /= norm; q02 /= norm; q03 /= norm;
    
    vector<double3> newTetramer2 = rotate(double4(q00,q01,q02,q03),refTetramer2);
    for (auto& t : newTetramer2) t = t + d1;
    
    const int iMax = ceil(rCut/a1l) + 1;
    const int jMax = ceil(rCut/a2l) + 1;
    const int kMax = ceil(rCut/a3l) + 1;
    
    double energy = 0.0;
    int tetramercounted = 0;
    vector<double> myeng(4,0.);
    
    //for (int i=-iMax; i<=iMax; i++)
    for (int i=-iMax; i<=iMax; i++)
    {
        
        for (int j=-jMax; j<=jMax; j++)
        {
            
            for (int k=-kMax; k<=kMax; k++)
            {
                const double3 site_ijk = double3(i*a1.x + j*a2.x + k*a3.x, i*a1.y + j*a2.y + k*a3.y, i*a1.z + j*a2.z + k*a3.z);
                //if (tetramercounted > 0) break;
                for (unsigned int l=0; l<newTetramer.size(); l++) {
                    double3 ref = newTetramer[l];
                    double3 ref2 = newTetramer2[l];
                    for (unsigned int m=0; m<newTetramer.size(); m++) {
                        double3 ijk = newTetramer[m] + site_ijk;
                        
                        double3 sijk = ref - ijk;
                        
                        /*std::cout << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                        std::cout << ref.x << "\t" << ref.y << "\t" << ref.z << "\n";
                        std::cout << sijk.x << "\t" << sijk.y << "\t" << sijk.z << "\n\n";*/
                        
                        
                        double r2 = sijk.x*sijk.x + sijk.y*sijk.y + sijk.z*sijk.z;
                        
                        if (r2 < rCut2 && r2 > 0.)
                        {
                            //energy += epsilon*exp(-r2/sigma2);
                            double ratio2 = sigma2/r2;
                            double ratio6 = pow(ratio2, 3.0);
                            double prefac = 1+alpha;
                            energy += 0.25*4.0*epsilon*prefac*ratio6*(ratio6 - 1.0);
                            myeng[l] += 4.0*epsilon*ratio6*(ratio6 - 1.0);
                        }
                        
                        
                        double3 ijk2 = newTetramer2[m] + site_ijk ;
                        sijk = ref - ijk2;
                        r2 = sijk.x*sijk.x + sijk.y*sijk.y + sijk.z*sijk.z;
                        
                        if (r2 < rCut2 && r2 > 0.)
                        {
                            //energy += epsilon*exp(-r2/sigma2);
                            double ratio2 = sigma2/r2;
                            double ratio6 = pow(ratio2, 3.0);
                            double prefac = 1-alpha;
                            energy += 2*0.25*4.0*epsilon*prefac*ratio6*(ratio6 - 1.0); //twice the contribution
                            myeng[l] += 2*4.0*epsilon*ratio6*(ratio6 - 1.0); //twice the contribution
                        }
                        
                        sijk = ref2 - ijk2;
                        r2 = sijk.x*sijk.x + sijk.y*sijk.y + sijk.z*sijk.z;
                        
                        if (r2 < rCut2 && r2 > 0.)
                        {
                            //energy += epsilon*exp(-r2/sigma2);
                            double ratio2 = sigma2/r2;
                            double ratio6 = pow(ratio2, 3.0);
                            double prefac = 1+alpha;
                            energy += 0.25*4.0*epsilon*prefac*ratio6*(ratio6 - 1.0);
                            myeng[l] += 4.0*epsilon*ratio6*(ratio6 - 1.0);
                        }

                        
                    }
                }
                tetramercounted+=2;
            }
        }
    }
    
    for (auto& i: myeng) i /= 2; //pairwise gives you (*2), four atoms per tetramer (/4), results in /= 2
    energy *= 2./2.; //pair wise, 2 (divided by 2 tetramers)
    
    const double fitness = exp(-energy);
    
    /*cout<<"Energy: "<<energy<<endl;
     cout<<"Fitness: "<<fitness<<endl;
     cin.get();*/
    
    return fitness;
}

namespace Racemic{int main(int, const char* []);};
namespace Enantiopure{int main(int, const char* []);};


int main(int argc, const char * argv[]) {
    bool racemic=true;
    if (argc == 3){
        if (atoi(argv[2])) racemic = true;
        else racemic = false;
        
    }
    if (racemic){
        alpha = -0.5;
        Racemic::main(argc,argv);
    }
    else{
        alpha = +0.5;
        Enantiopure::main(argc,argv);
    }
    
}

namespace Deprecate{
    int main(int argc, const char* argv[]){
        if (argc == 2)
        {
            density = atof(argv[1]);
        }
        cout<<"Calculating optimum lattice for Tetramer Model at density: "<<density<<endl;
        
        
        //initial config
        refTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
        refTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
        refTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
        refTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
        
        refTetramer2 = refTetramer;
        refTetramer2[3] = double3( 0.52915050,  0.26457525, 1.32287625);
        double3 translate = double3(0.,0.,0.);
        for (auto &i : refTetramer2) i = i + translate;
        
        // create initial generation
        vector<individual> individuals;
        const unsigned int populationCount = 1000;
        const unsigned int bit = 31;
        
        for (unsigned int i=0; i<populationCount; i++)
        {
            vector<gene> initGenes;
            
            //tetramer1 lattice positions
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI));
            
            //tetramer1 orientations
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            
            //tetramer2 positions
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            
            
            //tetramer2 orientations
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            
            
            individuals.push_back(individual(initGenes));
        }
        
        generation g = generation(individuals);
        g.randomize();
        g.elitism = true;
        
        ofstream fitnessData("fitness.dat");
        ofstream solutionParameterData("solutionParameter.dat");
        solutionParameterData<<"#g\tx\ty\tphi\tpsi\ttheta\tq0\tq1\tq2\tq3"<<endl;
        ofstream solutionVectorData("solutionVector.dat");
        solutionVectorData<<"#g\ta1x\ta1y\ta1z\ta2x\ta2y\ta2z\ta3x\ta3y\ta3z\trx\try\trz\tangle"<<endl;
        
        // calculate fitness of all individuals
        for (unsigned int t=0; t<500; t++)
        {
            g.calcFitness(fitness);
            g.sortByFitness();
            fitnessData<<t<<"\t"<<g.individuals[0].fitness<<endl;
            
            // Display and save properties of this generation's best individual
            for (unsigned int i=0; i<1; i++)
            {
                double q0,q1,q2,q3;
                double dx, dy, dz;
                double q00, q01, q02,q03;
                
                q0=g.individuals[i].genes[5].getValue();
                q1=g.individuals[i].genes[6].getValue();
                q2=g.individuals[i].genes[7].getValue();
                q3=g.individuals[i].genes[8].getValue();
                
                double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
                if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
                
                dx = g.individuals[i].genes[9].getValue();;
                dy = g.individuals[i].genes[10].getValue();
                dz = g.individuals[i].genes[11].getValue();
                
                q00=g.individuals[i].genes[12].getValue();
                q01=g.individuals[i].genes[13].getValue();
                q02=g.individuals[i].genes[14].getValue();
                q03=g.individuals[i].genes[15].getValue();
                
                norm = sqrt(q00*q00 + q01*q01 + q02*q02 + q03*q03);
                if (norm > 0.001) q00 /= norm; q01 /= norm; q02 /= norm; q03 /= norm;
                
                
                
                cout<<"fitness: "<<g.individuals[i].fitness<<", rel. fitness: "<<g.individuals[i].relativeFitness;
                cout<<", x: "<<g.individuals[i].genes[0].getValue();
                cout<<", y: "<<g.individuals[i].genes[1].getValue();
                cout<<", phi: "<<g.individuals[i].genes[2].getValue();
                cout<<", psi: "<<g.individuals[i].genes[3].getValue();
                cout<<", theta: "<<g.individuals[i].genes[4].getValue()<<endl;
                cout<<", q0: "<<q0<<endl;
                cout<<", q1: "<<q1<<endl;
                cout<<", q2: "<<q2<<endl;
                cout<<", q3: "<<q3<<endl;
                cout<<", dx: "<<dx<<endl;
                cout<<", dy: "<<dy<<endl;
                cout<<", dz: "<<dz<<endl;
                cout<<", q00: "<<q00<<endl;
                cout<<", q01: "<<q01<<endl;
                cout<<", q02: "<<q02<<endl;
                cout<<", q03: "<<q03<<endl;
                
                
                lattice ltc = lattice(g.individuals[i].genes[0].getValue(), g.individuals[i].genes[1].getValue(),
                                      g.individuals[i].genes[2].getValue(), g.individuals[i].genes[3].getValue(), g.individuals[i].genes[4].getValue());
                
                
                ltc.setDensity(density/2.);
                const double3 a1 = ltc.getLatticeVector(0);
                const double3 a2 = ltc.getLatticeVector(1);
                const double3 a3 = ltc.getLatticeVector(2);
                const double4 o1 = axisangle(double4(q0,q1,q2,q3));
                const double4 o2 = axisangle(double4(q00,q01,q02,q03));
                const double3 d1 = double3(a1.x*dx, a2.y*dy, a3.z*dz);
                
                
                cout<<"("<<a1.x<<", "<<a1.y<<", "<<a1.z<<")"<<endl;
                cout<<"("<<a2.x<<", "<<a2.y<<", "<<a2.z<<")"<<endl;
                cout<<"("<<a3.x<<", "<<a3.y<<", "<<a3.z<<")"<<endl;
                cout<<"("<<d1.x<<", "<<d1.y<<", "<<d1.z<<")"<<endl;
                cout<<"("<<o1.coords.x<<", "<<o1.coords.y<<", "<<o1.coords.z<<", "<<o1.scalar<<")"<<endl;
                cout<<"("<<o2.coords.x<<", "<<o2.coords.y<<", "<<o2.coords.z<<", "<<o2.scalar<<")"<<endl;
                
                
                // Save best solution to file
                solutionParameterData<<t<<"\t"<<g.individuals[i].genes[0].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[1].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[2].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[3].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[4].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[5].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[6].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[7].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[8].getValue()<<"\t";
                for (unsigned int j=9; j<15; j++) solutionParameterData<<g.individuals[i].genes[j].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[15].getValue()<<endl;
                
                solutionVectorData<<t<<"\t"<<a1.x<<"\t"<<a1.y<<"\t"<<a1.z<<"\t";
                solutionVectorData<<a2.x<<"\t"<<a2.y<<"\t"<<a2.z<<"\t";
                solutionVectorData<<a3.x<<"\t"<<a3.y<<"\t"<<a3.z<<"\t";
                solutionVectorData<<d1.x<<"\t"<<d1.y<<"\t"<<d1.z<<"\t";
                solutionVectorData<<o1.coords.x<<"\t"<<o1.coords.y<<"\t"<<o1.coords.z<<"\t"<<o1.scalar<<"\t";
                solutionVectorData<<o2.coords.x<<"\t"<<o2.coords.y<<"\t"<<o2.coords.z<<"\t"<<o2.scalar<<endl;
            }
            
            // create next generation
            g.mateFitnessProportional();
            g.mutate(1.0, 0.05);
            
            cout<<endl;
        }
        
        fitnessData.close();
        solutionParameterData.close();
        solutionVectorData.close();
        
        g.calcFitness(fitness);
        g.sortByFitness();
        std::cout << "fitness is: " << fitness(g.individuals[0]) << endl;
        writexyz(g.individuals[0]);
        
        
        return 0;
    }

}

namespace Racemic {
    void writexyz(individual _ind){
        // Create a lattice from the supplied individual
        lattice ltc = lattice(_ind.genes[0].getValue(), _ind.genes[1].getValue(), _ind.genes[2].getValue(), _ind.genes[3].getValue(), _ind.genes[4].getValue());
        ltc.setDensity(density/2.);
        
        double dx = _ind.genes[9].getValue();;
        double dy = _ind.genes[10].getValue();
        double dz = _ind.genes[11].getValue();
        
        
        const double3 a1 = ltc.getLatticeVector(0);
        const double3 a2 = ltc.getLatticeVector(1);
        const double3 a3 = ltc.getLatticeVector(2);
        const double a1l = sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
        const double a2l = sqrt(a2.x*a2.x + a2.y*a2.y + a2.z*a2.z);
        const double a3l = sqrt(a3.x*a3.x + a3.y*a3.y + a3.z*a3.z);
        const double3 d1 = double3(a1.x*dx, a2.y*dy, a3.z*dz);
        
        
        double q0,q1,q2,q3;
        q0=_ind.genes[5].getValue();
        q1=_ind.genes[6].getValue();
        q2=_ind.genes[7].getValue();
        q3=_ind.genes[8].getValue();
        
        double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
        if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
        
        vector<double3> newTetramer = rotate(double4(q0,q1,q2,q3),refTetramer);
        
        double q00, q01, q02,q03;
        q00=_ind.genes[12].getValue();
        q01=_ind.genes[13].getValue();
        q02=_ind.genes[14].getValue();
        q03=_ind.genes[15].getValue();
        
        norm = sqrt(q00*q00 + q01*q01 + q02*q02 + q03*q03);
        if (norm > 0.001) q00 /= norm; q01 /= norm; q02 /= norm; q03 /= norm;
        
        vector<double3> newTetramer2 = rotate(double4(q00,q01,q02,q03),refTetramer2);
        for (auto& t : newTetramer2) t = t + d1;
        
        const int iMax = ceil(rCut/a1l) + 1;
        const int jMax = ceil(rCut/a2l) + 1;
        const int kMax = ceil(rCut/a3l) + 1;
        
        ofstream myxyz("lattice.xyz");
        int tetramercounted = 0;
        int ncount=0;
        
        cout << "Box info\n";
        cout << 2.*rCut/a1l << "\t" << 2.*rCut/a2l << "\t" << 2.*rCut/a3l << "\n";
        cout << 2*iMax + 1 << "\t" << 2*jMax + 1 << "\t" << 2*kMax + 1 << "\n";
        myxyz << (2*iMax+1)*(2*jMax+1)*(2*kMax+1)*4 <<"\n\n";
        
        for (int i=-iMax; i<=iMax; i++)
        {
            for (int j=-jMax; j<=jMax; j++)
            {
                for (int k=-kMax; k<=kMax; k++)
                {
                    const double3 site_ijk = double3(i*a1.x + j*a2.x + k*a3.x, i*a1.y + j*a2.y + k*a3.y, i*a1.z + j*a2.z + k*a3.z);
                    if ((i==j) &&(j==k) && (k==0)){
                        cout << "I am up to " << ncount << "\n";
                    }
                    for (unsigned int l=0; l<newTetramer.size(); l++) {
                        double3 ijk = newTetramer[l] + site_ijk;
                        myxyz << "O\t" << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                        ncount +=3;
                    }
                    
                    for (unsigned int l=0; l<newTetramer2.size(); l++) {
                        double3 ijk = newTetramer2[l] + site_ijk; // + d1;
                        myxyz << "N\t" << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                        ncount +=3;
                    }
                    tetramercounted+=2;
                    
                }
            }
        }
        myxyz.close();
        ofstream myinfo("box_config.dat");
        //myinfo << "Box info for python script\n";
        myinfo << (2*iMax+1)*a1.x << "\t" << (2*jMax+1)*a2.x << "\t" << (2*jMax+1)*a2.y << "\t"
        << (2*kMax+1)*a3.x << "\t" << (2*kMax+1)*a3.y << "\t" << (2*kMax+1)*a3.z << "\n";
        myinfo.close();
        
    }
    
    double fitness(individual _ind)
    {
        // Create a lattice from the supplied individual
        lattice ltc = lattice(_ind.genes[0].getValue(), _ind.genes[1].getValue(), _ind.genes[2].getValue(), _ind.genes[3].getValue(), _ind.genes[4].getValue());
        
        ltc.setDensity(density/2.);
        
        double dx = _ind.genes[9].getValue();;
        double dy = _ind.genes[10].getValue();
        double dz = _ind.genes[11].getValue();
        
        
        const double3 a1 = ltc.getLatticeVector(0);
        const double3 a2 = ltc.getLatticeVector(1);
        const double3 a3 = ltc.getLatticeVector(2);
        const double a1l = sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
        const double a2l = sqrt(a2.x*a2.x + a2.y*a2.y + a2.z*a2.z);
        const double a3l = sqrt(a3.x*a3.x + a3.y*a3.y + a3.z*a3.z);
        const double3 d1 = double3(a1.x*dx, a2.y*dy, a3.z*dz);
        
        double q0,q1,q2,q3;
        q0=_ind.genes[5].getValue();
        q1=_ind.genes[6].getValue();
        q2=_ind.genes[7].getValue();
        q3=_ind.genes[8].getValue();
        
        double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
        if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
        
        vector<double3> newTetramer = rotate(double4(q0,q1,q2,q3),refTetramer);
        
        double q00, q01, q02,q03;
        q00=_ind.genes[12].getValue();
        q01=_ind.genes[13].getValue();
        q02=_ind.genes[14].getValue();
        q03=_ind.genes[15].getValue();
        
        norm = sqrt(q00*q00 + q01*q01 + q02*q02 + q03*q03);
        if (norm > 0.001) q00 /= norm; q01 /= norm; q02 /= norm; q03 /= norm;
        
        vector<double3> newTetramer2 = rotate(double4(q00,q01,q02,q03),refTetramer2);
        for (auto& t : newTetramer2) t = t + d1;
        
        const int iMax = ceil(rCut/a1l) + 1;
        const int jMax = ceil(rCut/a2l) + 1;
        const int kMax = ceil(rCut/a3l) + 1;
        
        double energy = 0.0;
        int tetramercounted = 0;
        vector<double> myeng(4,0.);
        
        //for (int i=-iMax; i<=iMax; i++)
        for (int i=-iMax; i<=iMax; i++)
        {
            
            for (int j=-jMax; j<=jMax; j++)
            {
                
                for (int k=-kMax; k<=kMax; k++)
                {
                    const double3 site_ijk = double3(i*a1.x + j*a2.x + k*a3.x, i*a1.y + j*a2.y + k*a3.y, i*a1.z + j*a2.z + k*a3.z);
                    //if (tetramercounted > 0) break;
                    for (unsigned int l=0; l<newTetramer.size(); l++) {
                        double3 ref = newTetramer[l];
                        double3 ref2 = newTetramer2[l];
                        for (unsigned int m=0; m<newTetramer.size(); m++) {
                            double3 ijk = newTetramer[m] + site_ijk;
                            
                            double3 sijk = ref - ijk;
                            
                            /*std::cout << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                             std::cout << ref.x << "\t" << ref.y << "\t" << ref.z << "\n";
                             std::cout << sijk.x << "\t" << sijk.y << "\t" << sijk.z << "\n\n";*/
                            
                            
                            double r2 = sijk.x*sijk.x + sijk.y*sijk.y + sijk.z*sijk.z;
                            
                            if (r2 < rCut2 && r2 > 0.)
                            {
                                //energy += epsilon*exp(-r2/sigma2);
                                double ratio2 = sigma2/r2;
                                double ratio6 = pow(ratio2, 3.0);
                                double prefac = 1+alpha;
                                energy += 0.25*4.0*epsilon*prefac*ratio6*(ratio6 - 1.0);
                                myeng[l] += 4.0*epsilon*ratio6*(ratio6 - 1.0);
                            }
                            
                            
                            double3 ijk2 = newTetramer2[m] + site_ijk ;
                            sijk = ref - ijk2;
                            r2 = sijk.x*sijk.x + sijk.y*sijk.y + sijk.z*sijk.z;
                            
                            if (r2 < rCut2 && r2 > 0.)
                            {
                                //energy += epsilon*exp(-r2/sigma2);
                                double ratio2 = sigma2/r2;
                                double ratio6 = pow(ratio2, 3.0);
                                double prefac = 1-alpha;
                                energy += 2*0.25*4.0*epsilon*prefac*ratio6*(ratio6 - 1.0); //twice the contribution
                                myeng[l] += 2*4.0*epsilon*ratio6*(ratio6 - 1.0); //twice the contribution
                            }
                            
                            sijk = ref2 - ijk2;
                            r2 = sijk.x*sijk.x + sijk.y*sijk.y + sijk.z*sijk.z;
                            
                            if (r2 < rCut2 && r2 > 0.)
                            {
                                //energy += epsilon*exp(-r2/sigma2);
                                double ratio2 = sigma2/r2;
                                double ratio6 = pow(ratio2, 3.0);
                                double prefac = 1+alpha;
                                energy += 0.25*4.0*epsilon*prefac*ratio6*(ratio6 - 1.0);
                                myeng[l] += 4.0*epsilon*ratio6*(ratio6 - 1.0);
                            }
                            
                            
                        }
                    }
                    tetramercounted+=2;
                }
            }
        }
        
        for (auto& i: myeng) i /= 2; //pairwise gives you (*2), four atoms per tetramer (/4), results in /= 2
        energy *= 2./2.; //pair wise, 2 (divided by 2 tetramers)
        
        const double fitness = exp(-energy);
        
        /*cout<<"Energy: "<<energy<<endl;
         cout<<"Fitness: "<<fitness<<endl;
         cin.get();*/
        
        return fitness;
    }
    
    int main(int argc, const char * argv[]) {
        if (argc >= 2)
        {
            density = atof(argv[1]);
        }
        cout<<"Calculating optimum lattice for Tetramer Model at density: "<<density<<endl;
        
        
        //initial config
        refTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
        refTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
        refTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
        refTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
        
        refTetramer2 = refTetramer;
        refTetramer2[3] = double3( 0.52915050,  0.26457525, 1.32287625);
        double3 translate = double3(0.,0.,0.);
        for (auto &i : refTetramer2) i = i + translate;
        
        // create initial generation
        vector<individual> individuals;
        const unsigned int populationCount = 1000;
        const unsigned int bit = 31;
        
        for (unsigned int i=0; i<populationCount; i++)
        {
            vector<gene> initGenes;
            
            //tetramer1 lattice positions
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI));
            
            //tetramer1 orientations
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            
            //tetramer2 positions
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            
            
            //tetramer2 orientations
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            
            
            individuals.push_back(individual(initGenes));
        }
        
        generation g = generation(individuals);
        g.randomize();
        g.elitism = true;
        
        ofstream fitnessData("fitness.dat");
        ofstream solutionParameterData("solutionParameter.dat");
        solutionParameterData<<"#g\tx\ty\tphi\tpsi\ttheta\tq0\tq1\tq2\tq3"<<endl;
        ofstream solutionVectorData("solutionVector.dat");
        solutionVectorData<<"#g\ta1x\ta1y\ta1z\ta2x\ta2y\ta2z\ta3x\ta3y\ta3z\trx\try\trz\tangle"<<endl;
        
        // calculate fitness of all individuals
        for (unsigned int t=0; t<500; t++)
        {
            g.calcFitness(fitness);
            g.sortByFitness();
            fitnessData<<t<<"\t"<<g.individuals[0].fitness<<endl;
            
            // Display and save properties of this generation's best individual
            for (unsigned int i=0; i<1; i++)
            {
                double q0,q1,q2,q3;
                double dx, dy, dz;
                double q00, q01, q02,q03;
                
                q0=g.individuals[i].genes[5].getValue();
                q1=g.individuals[i].genes[6].getValue();
                q2=g.individuals[i].genes[7].getValue();
                q3=g.individuals[i].genes[8].getValue();
                
                double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
                if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
                
                dx = g.individuals[i].genes[9].getValue();;
                dy = g.individuals[i].genes[10].getValue();
                dz = g.individuals[i].genes[11].getValue();
                
                q00=g.individuals[i].genes[12].getValue();
                q01=g.individuals[i].genes[13].getValue();
                q02=g.individuals[i].genes[14].getValue();
                q03=g.individuals[i].genes[15].getValue();
                
                norm = sqrt(q00*q00 + q01*q01 + q02*q02 + q03*q03);
                if (norm > 0.001) q00 /= norm; q01 /= norm; q02 /= norm; q03 /= norm;
                
                
                
                cout<<"fitness: "<<g.individuals[i].fitness<<", rel. fitness: "<<g.individuals[i].relativeFitness;
                cout<<", x: "<<g.individuals[i].genes[0].getValue();
                cout<<", y: "<<g.individuals[i].genes[1].getValue();
                cout<<", phi: "<<g.individuals[i].genes[2].getValue();
                cout<<", psi: "<<g.individuals[i].genes[3].getValue();
                cout<<", theta: "<<g.individuals[i].genes[4].getValue()<<endl;
                cout<<", q0: "<<q0<<endl;
                cout<<", q1: "<<q1<<endl;
                cout<<", q2: "<<q2<<endl;
                cout<<", q3: "<<q3<<endl;
                cout<<", dx: "<<dx<<endl;
                cout<<", dy: "<<dy<<endl;
                cout<<", dz: "<<dz<<endl;
                cout<<", q00: "<<q00<<endl;
                cout<<", q01: "<<q01<<endl;
                cout<<", q02: "<<q02<<endl;
                cout<<", q03: "<<q03<<endl;
                
                
                lattice ltc = lattice(g.individuals[i].genes[0].getValue(), g.individuals[i].genes[1].getValue(),
                                      g.individuals[i].genes[2].getValue(), g.individuals[i].genes[3].getValue(), g.individuals[i].genes[4].getValue());
                
                
                ltc.setDensity(density/2.);
                const double3 a1 = ltc.getLatticeVector(0);
                const double3 a2 = ltc.getLatticeVector(1);
                const double3 a3 = ltc.getLatticeVector(2);
                const double4 o1 = axisangle(double4(q0,q1,q2,q3));
                const double4 o2 = axisangle(double4(q00,q01,q02,q03));
                const double3 d1 = double3(a1.x*dx, a2.y*dy, a3.z*dz);
                
                
                cout<<"("<<a1.x<<", "<<a1.y<<", "<<a1.z<<")"<<endl;
                cout<<"("<<a2.x<<", "<<a2.y<<", "<<a2.z<<")"<<endl;
                cout<<"("<<a3.x<<", "<<a3.y<<", "<<a3.z<<")"<<endl;
                cout<<"("<<d1.x<<", "<<d1.y<<", "<<d1.z<<")"<<endl;
                cout<<"("<<o1.coords.x<<", "<<o1.coords.y<<", "<<o1.coords.z<<", "<<o1.scalar<<")"<<endl;
                cout<<"("<<o2.coords.x<<", "<<o2.coords.y<<", "<<o2.coords.z<<", "<<o2.scalar<<")"<<endl;
                
                
                // Save best solution to file
                solutionParameterData<<t<<"\t"<<g.individuals[i].genes[0].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[1].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[2].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[3].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[4].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[5].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[6].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[7].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[8].getValue()<<"\t";
                for (unsigned int j=9; j<15; j++) solutionParameterData<<g.individuals[i].genes[j].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[15].getValue()<<endl;
                
                solutionVectorData<<t<<"\t"<<a1.x<<"\t"<<a1.y<<"\t"<<a1.z<<"\t";
                solutionVectorData<<a2.x<<"\t"<<a2.y<<"\t"<<a2.z<<"\t";
                solutionVectorData<<a3.x<<"\t"<<a3.y<<"\t"<<a3.z<<"\t";
                solutionVectorData<<d1.x<<"\t"<<d1.y<<"\t"<<d1.z<<"\t";
                solutionVectorData<<o1.coords.x<<"\t"<<o1.coords.y<<"\t"<<o1.coords.z<<"\t"<<o1.scalar<<"\t";
                solutionVectorData<<o2.coords.x<<"\t"<<o2.coords.y<<"\t"<<o2.coords.z<<"\t"<<o2.scalar<<endl;
            }
            
            // create next generation
            g.mateFitnessProportional();
            g.mutate(1.0, 0.05);
            
            cout<<endl;
        }
        
        fitnessData.close();
        solutionParameterData.close();
        solutionVectorData.close();
        
        g.calcFitness(Racemic::fitness);
        g.sortByFitness();
        std::cout << "fitness is: " << Racemic::fitness(g.individuals[0]) << endl;
        Racemic::writexyz(g.individuals[0]);
        
        
        return 0;
        
    }
}


namespace Enantiopure{
    
    void writexyz(individual _ind){
        // Create a lattice from the supplied individual
        lattice ltc = lattice(_ind.genes[0].getValue(), _ind.genes[1].getValue(), _ind.genes[2].getValue(), _ind.genes[3].getValue(), _ind.genes[4].getValue());
        ltc.setDensity(density);
        
        const double3 a1 = ltc.getLatticeVector(0);
        const double3 a2 = ltc.getLatticeVector(1);
        const double3 a3 = ltc.getLatticeVector(2);
        const double a1l = sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
        const double a2l = sqrt(a2.x*a2.x + a2.y*a2.y + a2.z*a2.z);
        const double a3l = sqrt(a3.x*a3.x + a3.y*a3.y + a3.z*a3.z);
        
        double q0,q1,q2,q3;
        q0=_ind.genes[5].getValue();
        q1=_ind.genes[6].getValue();
        q2=_ind.genes[7].getValue();
        q3=_ind.genes[8].getValue();
        
        double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
        if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
        
        vector<double3> newTetramer = rotate(double4(q0,q1,q2,q3),refTetramer);
    
        
        const int iMax = ceil(rCut/a1l) + 1;
        const int jMax = ceil(rCut/a2l) + 1;
        const int kMax = ceil(rCut/a3l) + 1;
        
        ofstream myxyz("lattice.xyz");
        int tetramercounted = 0;
        int ncount=0;
        
        cout << "Box info\n";
        cout << 2.*rCut/a1l << "\t" << 2.*rCut/a2l << "\t" << 2.*rCut/a3l << "\n";
        cout << 2*iMax + 1 << "\t" << 2*jMax + 1 << "\t" << 2*kMax + 1 << "\n";
        myxyz << (2*iMax+1)*(2*jMax+1)*(2*kMax+1)*4 <<"\n\n";
        
        for (int i=-iMax; i<=iMax; i++)
        {
            for (int j=-jMax; j<=jMax; j++)
            {
                for (int k=-kMax; k<=kMax; k++)
                {
                    const double3 site_ijk = double3(i*a1.x + j*a2.x + k*a3.x, i*a1.y + j*a2.y + k*a3.y, i*a1.z + j*a2.z + k*a3.z);
                    if ((i==j) &&(j==k) && (k==0)){
                        cout << "I am up to " << ncount << "\n";
                    }
                    for (unsigned int l=0; l<newTetramer.size(); l++) {
                        double3 ijk = newTetramer[l] + site_ijk;
                        myxyz << "O\t" << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                        ncount +=3;
                    }
                    tetramercounted++;
                    
                }
            }
        }
        myxyz.close();
        ofstream myinfo("box_config.dat");
        myinfo << "Box info for python script\n";
        myinfo << (2*iMax+1)*a1.x << "\t" << (2*jMax+1)*a2.x << "\t" << (2*jMax+1)*a2.y << "\t"
             << (2*kMax+1)*a3.x << "\t" << (2*kMax+1)*a3.y << "\t" << (2*kMax+1)*a3.z << "\n";
        myinfo.close();
        
    }

    
    double fitness(individual _ind)
    {
        // Create a lattice from the supplied individual
        lattice ltc = lattice(_ind.genes[0].getValue(), _ind.genes[1].getValue(), _ind.genes[2].getValue(), _ind.genes[3].getValue(), _ind.genes[4].getValue());
        ltc.setDensity(density);
        
        const double3 a1 = ltc.getLatticeVector(0);
        const double3 a2 = ltc.getLatticeVector(1);
        const double3 a3 = ltc.getLatticeVector(2);
        const double a1l = sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
        const double a2l = sqrt(a2.x*a2.x + a2.y*a2.y + a2.z*a2.z);
        const double a3l = sqrt(a3.x*a3.x + a3.y*a3.y + a3.z*a3.z);
        
        double q0,q1,q2,q3;
        q0=_ind.genes[5].getValue();
        q1=_ind.genes[6].getValue();
        q2=_ind.genes[7].getValue();
        q3=_ind.genes[8].getValue();
        
        double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
        if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
        
        vector<double3> newTetramer = rotate(double4(q0,q1,q2,q3),refTetramer);
        
        const int iMax = ceil(rCut/a1l) + 1;
        const int jMax = ceil(rCut/a2l) + 1;
        const int kMax = ceil(rCut/a3l) + 1;
        
        double energy = 0.0;
        int tetramercounted = 0;
        vector<double> myeng(4,0.);
        
        //for (int i=-iMax; i<=iMax; i++)
        for (int i=-iMax; i<=iMax; i++)
        {
            
            for (int j=-jMax; j<=jMax; j++)
            {
                for (int k=-kMax; k<=kMax; k++)
                {
                    const double3 site_ijk = double3(i*a1.x + j*a2.x + k*a3.x, i*a1.y + j*a2.y + k*a3.y, i*a1.z + j*a2.z + k*a3.z);
                    //if (tetramercounted > 0) break;
                    for (unsigned int l=0; l<newTetramer.size(); l++) {
                        double3 ref = newTetramer[l];
                        for (unsigned int m=0; m<newTetramer.size(); m++) {
                            double3 ijk = newTetramer[m] + site_ijk;
                            
                            double3 sijk = ref - ijk;
                            
                            /*std::cout << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                             std::cout << ref.x << "\t" << ref.y << "\t" << ref.z << "\n";
                             std::cout << sijk.x << "\t" << sijk.y << "\t" << sijk.z << "\n\n";*/
                            
                            
                            const double r2 = sijk.x*sijk.x + sijk.y*sijk.y + sijk.z*sijk.z;
                            
                            if (r2 < rCut2 && r2 > 0.)
                            {
                                //energy += epsilon*exp(-r2/sigma2);
                                double ratio2 = sigma2/r2;
                                double ratio6 = pow(ratio2, 3.0);
                                energy += 0.25*4.0*(1+alpha)*epsilon*ratio6*(ratio6 - 1.0);
                                myeng[l] += 4.0*epsilon*(1+alpha)*ratio6*(ratio6 - 1.0);
                            }
                        }
                    }
                    tetramercounted++;
                }
            }
        }
        
        for (auto& i: myeng) i /= 2; //pairwise gives you (*2), four atoms per tetramer (/4), results in /= 2
        energy *= 2.; //pair wise, 2
        
        const double fitness = exp(-energy);
        
        /*cout<<"Energy: "<<energy<<endl;
         cout<<"Fitness: "<<fitness<<endl;
         cin.get();*/
        
        return fitness;
    }

    
    int main(int argc, const char * argv[]) {
        if (argc >= 2)
        {
            density = atof(argv[1]);
        }
        cout<<"Calculating optimum lattice for Tetramer Model at density: "<<density<<endl;
        
        
        //initial config
        refTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
        refTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
        refTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
        refTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
        
        // create initial generation
        vector<individual> individuals;
        const unsigned int populationCount = 1000;
        const unsigned int bit = 31;
        
        for (unsigned int i=0; i<populationCount; i++)
        {
            vector<gene> initGenes;
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI));
            
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            
            
            
            individuals.push_back(individual(initGenes));
        }
        
        generation g = generation(individuals);
        g.randomize();
        g.elitism = true;
        
        ofstream fitnessData("fitness.dat");
        ofstream solutionParameterData("solutionParameter.dat");
        solutionParameterData<<"#g\tx\ty\tphi\tpsi\ttheta\tq0\tq1\tq2\tq3"<<endl;
        ofstream solutionVectorData("solutionVector.dat");
        solutionVectorData<<"#g\ta1x\ta1y\ta1z\ta2x\ta2y\ta2z\ta3x\ta3y\ta3z\trx\try\trz\tangle"<<endl;
        
        // calculate fitness of all individuals
        for (unsigned int t=0; t<1000; t++)
        {
            g.calcFitness(Enantiopure::fitness);
            g.sortByFitness();
            fitnessData<<t<<"\t"<<g.individuals[0].fitness<<endl;
            
            // Display and save properties of this generation's best individual
            for (unsigned int i=0; i<1; i++)
            {
                double q0,q1,q2,q3;
                q0=g.individuals[i].genes[5].getValue();
                q1=g.individuals[i].genes[6].getValue();
                q2=g.individuals[i].genes[7].getValue();
                q3=g.individuals[i].genes[8].getValue();
                
                
                double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
                if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
                
                
                cout<<"fitness: "<<g.individuals[i].fitness<<", rel. fitness: "<<g.individuals[i].relativeFitness;
                cout<<", x: "<<g.individuals[i].genes[0].getValue();
                cout<<", y: "<<g.individuals[i].genes[1].getValue();
                cout<<", phi: "<<g.individuals[i].genes[2].getValue();
                cout<<", psi: "<<g.individuals[i].genes[3].getValue();
                cout<<", theta: "<<g.individuals[i].genes[4].getValue()<<endl;
                cout<<", q0:"<<q0<<endl;
                cout<<", q1:"<<q1<<endl;
                cout<<", q2:"<<q2<<endl;
                cout<<", q3:"<<q3<<endl;
                
                
                lattice ltc = lattice(g.individuals[i].genes[0].getValue(), g.individuals[i].genes[1].getValue(),
                                      g.individuals[i].genes[2].getValue(), g.individuals[i].genes[3].getValue(), g.individuals[i].genes[4].getValue());
                
                ltc.setDensity(density);
                const double3 a1 = ltc.getLatticeVector(0);
                const double3 a2 = ltc.getLatticeVector(1);
                const double3 a3 = ltc.getLatticeVector(2);
                const double4 o1 = axisangle(double4(q0,q1,q2,q3));
                
                
                cout<<"("<<a1.x<<", "<<a1.y<<", "<<a1.z<<")"<<endl;
                cout<<"("<<a2.x<<", "<<a2.y<<", "<<a2.z<<")"<<endl;
                cout<<"("<<a3.x<<", "<<a3.y<<", "<<a3.z<<")"<<endl;
                cout<<"("<<o1.coords.x<<", "<<o1.coords.y<<", "<<o1.coords.z<<", "<<o1.scalar<<")"<<endl;
                
                
                // Save best solution to file
                solutionParameterData<<t<<"\t"<<g.individuals[i].genes[0].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[1].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[2].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[3].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[4].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[5].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[6].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[7].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[8].getValue()<<endl;
                
                solutionVectorData<<t<<"\t"<<a1.x<<"\t"<<a1.y<<"\t"<<a1.z<<"\t";
                solutionVectorData<<a2.x<<"\t"<<a2.y<<"\t"<<a2.z<<"\t";
                solutionVectorData<<a3.x<<"\t"<<a3.y<<"\t"<<a3.z<<"\t";
                solutionVectorData<<o1.coords.x<<"\t"<<o1.coords.y<<"\t"<<o1.coords.z<<"\t"<<o1.scalar<<endl;
            }
            
            // create next generation
            g.mateFitnessProportional();
            g.mutate(1.0, 0.05);
            
            cout<<endl;
        }
        
        fitnessData.close();
        solutionParameterData.close();
        solutionVectorData.close();
        
        g.calcFitness(Enantiopure::fitness);
        g.sortByFitness();
        std::cout << "fitness is: " << Enantiopure::fitness(g.individuals[0]) << endl;
        Enantiopure::writexyz(g.individuals[0]);
        
        
        return 0;

    }
    
    
}

namespace old2 {
    int main(int argc, char* argv[])
    {
        double mycount=2;
        density=0.24;
        vector<double> _rvec(9,0.);
        _rvec[0]= 0.599691; _rvec[1] = 0.793268; _rvec[2] = 1.2481; _rvec[3] = 1.56853; _rvec[4]= 1.57546;
        _rvec[5]= -0.338845; _rvec[6] = 0.313275; _rvec[7] = -0.888193; _rvec[8] = -1.27391e-05;
        
        
        //initail config
        refTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
        refTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
        refTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
        refTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
    
        
        
        
        
        // create initial generation
        vector<individual> individuals;
        const unsigned int populationCount = mycount;
        const unsigned int bit = 32;
        
        for (unsigned int i=0; i<populationCount; i++)
        {
            vector<gene> initGenes;
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI));
            
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            
            
            
            individuals.push_back(individual(initGenes));
        }
        
        generation g = generation(individuals);
        g.randomize();
        g.elitism = true;
        
        for (unsigned int ii=0; ii<_rvec.size(); ii++) g.individuals[0].genes[ii].setValue(_rvec[ii]);
        
        ofstream fitnessData("fitness.dat");
        ofstream solutionParameterData("solutionParameter.dat");
        solutionParameterData<<"#g\tx\ty\tphi\tpsi\ttheta\tq0\tq1\tq2\tq3"<<endl;
        ofstream solutionVectorData("solutionVector.dat");
        solutionVectorData<<"#g\ta1x\ta1y\ta1z\ta2x\ta2y\ta2z\ta3x\ta3y\ta3z\trx\try\trz\tangle"<<endl;
        
        // calculate fitness of all individuals
        for (unsigned int t=0; t<mycount; t++)
        {
            g.calcFitness(fitness);
            g.sortByFitness();
            fitnessData<<t<<"\t"<<g.individuals[0].fitness<<endl;
            
            // Display and save properties of this generation's best individual
            for (unsigned int i=0; i<1; i++)
            {
                double q0,q1,q2,q3;
                q0=g.individuals[i].genes[5].getValue();
                q1=g.individuals[i].genes[6].getValue();
                q2=g.individuals[i].genes[7].getValue();
                q3=g.individuals[i].genes[8].getValue();
                
                
                double norm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
                if (norm > 0.001) q0 /= norm; q1 /= norm; q2 /= norm; q3 /= norm;
                
                
                cout<<"fitness: "<<g.individuals[i].fitness<<", rel. fitness: "<<g.individuals[i].relativeFitness;
                cout<<", x: "<<g.individuals[i].genes[0].getValue();
                cout<<", y: "<<g.individuals[i].genes[1].getValue();
                cout<<", phi: "<<g.individuals[i].genes[2].getValue();
                cout<<", psi: "<<g.individuals[i].genes[3].getValue();
                cout<<", theta: "<<g.individuals[i].genes[4].getValue()<<endl;
                cout<<", q0:"<<q0<<endl;
                cout<<", q1:"<<q1<<endl;
                cout<<", q2:"<<q2<<endl;
                cout<<", q3:"<<q3<<endl;
                
                
                lattice ltc = lattice(g.individuals[i].genes[0].getValue(), g.individuals[i].genes[1].getValue(),
                                      g.individuals[i].genes[2].getValue(), g.individuals[i].genes[3].getValue(), g.individuals[i].genes[4].getValue());
                
                ltc.setDensity(density);
                const double3 a1 = ltc.getLatticeVector(0);
                const double3 a2 = ltc.getLatticeVector(1);
                const double3 a3 = ltc.getLatticeVector(2);
                const double4 o1 = axisangle(double4(q0,q1,q2,q3));
                
                
                cout<<"("<<a1.x<<", "<<a1.y<<", "<<a1.z<<")"<<endl;
                cout<<"("<<a2.x<<", "<<a2.y<<", "<<a2.z<<")"<<endl;
                cout<<"("<<a3.x<<", "<<a3.y<<", "<<a3.z<<")"<<endl;
                cout<<"("<<o1.coords.x<<", "<<o1.coords.y<<", "<<o1.coords.z<<", "<<o1.scalar<<")"<<endl;
                
                
                // Save best solution to file
                solutionParameterData<<t<<"\t"<<g.individuals[i].genes[0].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[1].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[2].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[3].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[4].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[5].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[6].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[7].getValue()<<"\t";
                solutionParameterData<<g.individuals[i].genes[8].getValue()<<endl;
                
                solutionVectorData<<t<<"\t"<<a1.x<<"\t"<<a1.y<<"\t"<<a1.z<<"\t";
                solutionVectorData<<a2.x<<"\t"<<a2.y<<"\t"<<a2.z<<"\t";
                solutionVectorData<<a3.x<<"\t"<<a3.y<<"\t"<<a3.z<<"\t";
                solutionVectorData<<o1.coords.x<<"\t"<<o1.coords.y<<"\t"<<o1.coords.z<<"\t"<<o1.scalar<<endl;
            }
            
            // create next generation
            g.mateFitnessProportional();
            g.mutate(1.0, 0.05);
            
            cout<<endl;
        }
        
        fitnessData.close();
        solutionParameterData.close();
        solutionVectorData.close();
        
        g.calcFitness(fitness);
        g.sortByFitness();
        writexyz(g.individuals[0]);
        
        
        return 0;
    }


}

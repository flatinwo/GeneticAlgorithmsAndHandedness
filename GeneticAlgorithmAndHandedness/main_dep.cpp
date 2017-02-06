//
//  main_dep.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 12/29/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include <stdio.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "lattice.hpp"
#include "gene.hpp"
#include "individual.hpp"
#include "generation.hpp"
#include "operators.hpp"
#include "crystal_search.hpp"

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

//3.71375e+17

namespace Racemic {
    void writexyz(individual _ind){
        // Create a lattice from the supplied individual
        lattice ltc = lattice(_ind.genes[0].getValue(), _ind.genes[1].getValue(), _ind.genes[2].getValue(), _ind.genes[3].getValue(), _ind.genes[4].getValue());
        ltc.setDensity(density/4.);
        
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
        zeroCOM(newTetramer);
        
        double q00, q01, q02,q03;
        q00=_ind.genes[12].getValue();
        q01=_ind.genes[13].getValue();
        q02=_ind.genes[14].getValue();
        q03=_ind.genes[15].getValue();
        
        norm = sqrt(q00*q00 + q01*q01 + q02*q02 + q03*q03);
        if (norm > 0.001) q00 /= norm; q01 /= norm; q02 /= norm; q03 /= norm;
        vector<double3> newTetramer2 = rotate(double4(q00,q01,q02,q03),refTetramer2);
        zeroCOM(newTetramer2);
        for (auto& t : newTetramer2) t = t + d1;
        
        double dx1 = _ind.genes[16].getValue();;
        double dy1 = _ind.genes[17].getValue();
        double dz1 = _ind.genes[18].getValue();
        const double3 d2 = double3(a1.x*dx1, a2.y*dy1, a3.z*dz1);
        
        double q10=_ind.genes[19].getValue();
        double q11=_ind.genes[20].getValue();
        double q12=_ind.genes[21].getValue();
        double q13=_ind.genes[22].getValue();
        
        norm = sqrt(q10*q10 + q11*q11 + q12*q12 + q13*q13);
        if (norm > 0.001) q10 /= norm; q11 /= norm; q12 /= norm; q13 /= norm;
        vector<double3> newTetramer3 = rotate(double4(q10,q11,q12,q13),refTetramer);
        zeroCOM(newTetramer3);
        for (auto& t : newTetramer3) t = t + d2;
        
        double dx2 = _ind.genes[23].getValue();;
        double dy2 = _ind.genes[24].getValue();
        double dz2 = _ind.genes[25].getValue();
        const double3 d3 = double3(a1.x*dx2, a2.y*dy2, a3.z*dz2);
        
        double q20=_ind.genes[26].getValue();
        double q21=_ind.genes[27].getValue();
        double q22=_ind.genes[28].getValue();
        double q23=_ind.genes[29].getValue();
        
        norm = sqrt(q20*q20 + q21*q21 + q22*q22 + q23*q23);
        if (norm > 0.001) q20 /= norm; q21 /= norm; q22 /= norm; q23 /= norm;
        vector<double3> newTetramer4 = rotate(double4(q20,q21,q22,q23),refTetramer2);
        zeroCOM(newTetramer4);
        for (auto& t : newTetramer4) t = t + d3;
        
        
        const double myrCut = 4.0;
        const int iMax = ceil(myrCut/a1l) + 1;
        const int jMax = ceil(myrCut/a2l) + 1;
        const int kMax = ceil(myrCut/a3l) + 1;
        
        ofstream myxyz("lattice.xyz");
        int tetramercounted = 0;
        int ncount=0;
        
        cout << "Box info\n";
        cout << 2.*myrCut/a1l << "\t" << 2.*myrCut/a2l << "\t" << 2.*myrCut/a3l << "\n";
        cout << 2*iMax + 1 << "\t" << 2*jMax + 1 << "\t" << 2*kMax + 1 << "\n";
        myxyz << (2*iMax+1)*(2*jMax+1)*(2*kMax+1)*8 <<"\n\n";
        
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
                        myxyz << "O1\t" << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                        ncount +=3;
                    }
                    
                    for (unsigned int l=0; l<newTetramer2.size(); l++) {
                        double3 ijk = newTetramer2[l] + site_ijk; // + d1;
                        myxyz << "N1\t" << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                        ncount +=3;
                    }
                    
                    for (unsigned int l=0; l<newTetramer.size(); l++) {
                        double3 ijk = newTetramer3[l] + site_ijk;
                        myxyz << "O2\t" << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                        ncount +=3;
                    }
                    
                    for (unsigned int l=0; l<newTetramer.size(); l++) {
                        double3 ijk = newTetramer4[l] + site_ijk; // + d1;
                        myxyz << "N2\t" << ijk.x << "\t" << ijk.y << "\t" << ijk.z << "\n";
                        ncount +=3;
                    }
                    tetramercounted+=4;
                    
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
        
        ltc.setDensity(density/4.);
        
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
        zeroCOM(newTetramer);
        
        double q00, q01, q02,q03;
        q00=_ind.genes[12].getValue();
        q01=_ind.genes[13].getValue();
        q02=_ind.genes[14].getValue();
        q03=_ind.genes[15].getValue();
        
        norm = sqrt(q00*q00 + q01*q01 + q02*q02 + q03*q03);
        if (norm > 0.001) q00 /= norm; q01 /= norm; q02 /= norm; q03 /= norm;
        vector<double3> newTetramer2 = rotate(double4(q00,q01,q02,q03),refTetramer2);
        zeroCOM(newTetramer2);
        for (auto& t : newTetramer2) t = t + d1;
        
        double dx1 = _ind.genes[16].getValue();;
        double dy1 = _ind.genes[17].getValue();
        double dz1 = _ind.genes[18].getValue();
        const double3 d2 = double3(a1.x*dx1, a2.y*dy1, a3.z*dz1);
        
        double q10=_ind.genes[19].getValue();
        double q11=_ind.genes[20].getValue();
        double q12=_ind.genes[21].getValue();
        double q13=_ind.genes[22].getValue();
        
        norm = sqrt(q10*q10 + q11*q11 + q12*q12 + q13*q13);
        if (norm > 0.001) q10 /= norm; q11 /= norm; q12 /= norm; q13 /= norm;
        vector<double3> newTetramer3 = rotate(double4(q10,q11,q12,q13),refTetramer);
        zeroCOM(newTetramer3);
        for (auto& t : newTetramer3) t = t + d2;
        
        double dx2 = _ind.genes[23].getValue();;
        double dy2 = _ind.genes[24].getValue();
        double dz2 = _ind.genes[25].getValue();
        const double3 d3 = double3(a1.x*dx2, a2.y*dy2, a3.z*dz2);
        
        double q20=_ind.genes[26].getValue();
        double q21=_ind.genes[27].getValue();
        double q22=_ind.genes[28].getValue();
        double q23=_ind.genes[29].getValue();
        
        norm = sqrt(q20*q20 + q21*q21 + q22*q22 + q23*q23);
        if (norm > 0.001) q20 /= norm; q21 /= norm; q22 /= norm; q23 /= norm;
        vector<double3> newTetramer4 = rotate(double4(q20,q21,q22,q23),refTetramer2);
        zeroCOM(newTetramer4);
        for (auto& t : newTetramer4) t = t + d3;
        
        
        
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
                        double3 ref3 = newTetramer3[l];
                        double3 ref4 = newTetramer4[l];
                        for (unsigned int m=0; m<newTetramer.size(); m++) {
                            double3 ijk = newTetramer[m] + site_ijk;
                            
                            
                            //1-1
                            double3 sijk = ref - ijk;
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
                            
                            //1-2
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
                            
                            
                            //2-2
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
                            
                            
                            
                            //1-3
                            ijk2 = newTetramer3[m] + site_ijk ;
                            sijk = ref - ijk2;
                            r2 = sijk.x*sijk.x + sijk.y*sijk.y + sijk.z*sijk.z;
                            
                            if (r2 < rCut2 && r2 > 0.)
                            {
                                //energy += epsilon*exp(-r2/sigma2);
                                double ratio2 = sigma2/r2;
                                double ratio6 = pow(ratio2, 3.0);
                                double prefac = 1+alpha;
                                energy += 2*0.25*4.0*epsilon*prefac*ratio6*(ratio6 - 1.0); //twice the contribution
                                myeng[l] += 2*4.0*epsilon*ratio6*(ratio6 - 1.0); //twice the contribution
                            }
                            
                            //1-4
                            ijk2 = newTetramer4[m] + site_ijk ;
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
                            
                            
                            //2-3
                            ijk2 = newTetramer3[m] + site_ijk ;
                            sijk = ref2 - ijk2;
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
                            
                            //2-4
                            ijk2 = newTetramer4[m] + site_ijk ;
                            sijk = ref2 - ijk2;
                            r2 = sijk.x*sijk.x + sijk.y*sijk.y + sijk.z*sijk.z;
                            
                            if (r2 < rCut2 && r2 > 0.)
                            {
                                //energy += epsilon*exp(-r2/sigma2);
                                double ratio2 = sigma2/r2;
                                double ratio6 = pow(ratio2, 3.0);
                                double prefac = 1+alpha;
                                energy += 2*0.25*4.0*epsilon*prefac*ratio6*(ratio6 - 1.0); //twice the contribution
                                myeng[l] += 2*4.0*epsilon*ratio6*(ratio6 - 1.0); //twice the contribution
                            }
                            
                            //3-3
                            ijk2 = newTetramer3[m] + site_ijk ;
                            sijk = ref3 - ijk2;
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
                            
                            
                            //3-4
                            ijk2 = newTetramer4[m] + site_ijk ;
                            sijk = ref3 - ijk2;
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
                            
                            //4-4
                            ijk2 = newTetramer4[m] + site_ijk ;
                            sijk = ref4 - ijk2;
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
                    tetramercounted+=4;
                }
            }
        }
        
        for (auto& i: myeng) i /= 2; //pairwise gives you (*2), four atoms per tetramer (/4), results in /= 2
        energy *= 2./4.; //pair wise, 2 (divided by 2 tetramers)
        
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
        density=0.26;
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
            
            
            //tetramer3 positions
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            
            
            //tetramer3 orientations
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            initGenes.push_back(gene(bit,-1.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            
            
            //tetramer4 positions
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            
            
            //tetramer4 orientations
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
                
                
                double dx1 = g.individuals[i].genes[16].getValue();;
                double dy1 = g.individuals[i].genes[17].getValue();
                double dz1 = g.individuals[i].genes[18].getValue();
                
                double q10=g.individuals[i].genes[19].getValue();
                double q11=g.individuals[i].genes[20].getValue();
                double q12=g.individuals[i].genes[21].getValue();
                double q13=g.individuals[i].genes[22].getValue();
                
                norm = sqrt(q10*q10 + q11*q11 + q12*q12 + q13*q13);
                if (norm > 0.001) q10 /= norm; q11 /= norm; q12 /= norm; q13 /= norm;
                
                double dx2 = g.individuals[i].genes[23].getValue();;
                double dy2 = g.individuals[i].genes[24].getValue();
                double dz2 = g.individuals[i].genes[25].getValue();
                
                double q20=g.individuals[i].genes[26].getValue();
                double q21=g.individuals[i].genes[27].getValue();
                double q22=g.individuals[i].genes[28].getValue();
                double q23=g.individuals[i].genes[29].getValue();
                
                norm = sqrt(q20*q20 + q21*q21 + q22*q22 + q23*q23);
                if (norm > 0.001) q20 /= norm; q21 /= norm; q22 /= norm; q23 /= norm;
                
                
                
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
                cout<<", dx1: "<<dx1<<endl;
                cout<<", dy1: "<<dy1<<endl;
                cout<<", dz1: "<<dz1<<endl;
                cout<<", q10: "<<q10<<endl;
                cout<<", q11: "<<q11<<endl;
                cout<<", q12: "<<q12<<endl;
                cout<<", q13: "<<q13<<endl;
                cout<<", dx2: "<<dx2<<endl;
                cout<<", dy2: "<<dy2<<endl;
                cout<<", dz2: "<<dz2<<endl;
                cout<<", q20: "<<q20<<endl;
                cout<<", q21: "<<q21<<endl;
                cout<<", q22: "<<q22<<endl;
                cout<<", q23: "<<q23<<endl;
                
                
                lattice ltc = lattice(g.individuals[i].genes[0].getValue(), g.individuals[i].genes[1].getValue(),
                                      g.individuals[i].genes[2].getValue(), g.individuals[i].genes[3].getValue(), g.individuals[i].genes[4].getValue());
                
                
                ltc.setDensity(density/4.);
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
        q0=2.*_ind.genes[5].getValue()-1.;
        q1=2.*_ind.genes[6].getValue()-1.;
        q2=2.*_ind.genes[7].getValue()-1.;
        q3=2.*_ind.genes[8].getValue()-1.;
        
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
        q0=2.*_ind.genes[5].getValue()-1.;
        q1=2.*_ind.genes[6].getValue()-1.;
        q2=2.*_ind.genes[7].getValue()-1.;
        q3=2.*_ind.genes[8].getValue()-1.;
        
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
        const unsigned int bit = 32;
        
        for (unsigned int i=0; i<populationCount; i++)
        {
            vector<gene> initGenes;
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, 1.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI/2.0));
            initGenes.push_back(gene(bit, 0, PI));
            
            initGenes.push_back(gene(bit,0.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,0.0,1.0)); // for quaternions
            initGenes.push_back(gene(bit,0.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            initGenes.push_back(gene(bit,0.0,1.0)); // for quaternions, fourth is chosen to satisfy constraint
            
            
            
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
                q0=2*g.individuals[i].genes[5].getValue()-1.;
                q1=2*g.individuals[i].genes[6].getValue()-1.;
                q2=2*g.individuals[i].genes[7].getValue()-1.;
                q3=2*g.individuals[i].genes[8].getValue()-1.;
                
                
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

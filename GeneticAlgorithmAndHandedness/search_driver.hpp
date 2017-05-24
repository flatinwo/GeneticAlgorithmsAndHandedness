//
//  search_driver.hpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/23/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#ifndef search_driver_hpp
#define search_driver_hpp

#include <stdio.h>
#include "crystal_search.hpp"
#include <vector>
#include <random>

//encapsulate here for connecting with boost python

struct randomstore{
    
    unsigned min, max;
    randomstore(unsigned int a=2, unsigned int b=10):min(a),max(b){
    }
    
    unsigned generate(){
        return std::uniform_int_distribution<unsigned>{min,max}(eng);
    }    
    
    unsigned generate(unsigned _min, unsigned _max){
        return std::uniform_int_distribution<unsigned>{_min,_max}(eng);
    }
    
private:
    std::mt19937 eng{std::random_device{}()};
};

class SearchDriver{
public:
    
    enum SEARCHTYPE { ENANTIOPURE=0, RACEMIC=1, RACEMIZING=2};
    enum BITTYPE {RANDOMSWEEP=0, FIXEDSWEEP=1, NOSWEEP=2, NARROWINGRANDOMSWEEP=3};
    
    SearchDriver(short, SEARCHTYPE,
                 BITTYPE,
                 double density,
                 double lambda,
                 int popCount,
                 int iterations,
                 int bitmin,
                 int bitmax);
    
    ~SearchDriver();
    
    virtual void run();
    
    void setSeed(unsigned seed);
    void setBitType(BITTYPE);
    void updateBittage(unsigned int);
    void setIterationSweep(unsigned int);
 
    double getEnergy(){return _search->getLatticeEnergy();};
    
    void overrideGenes(std::vector<double>&);
    void overrideGenes(){ _search->overrideGeneration();};
    

protected:
    const int _nmol;
    std::vector<molecule> _molecule_list;
    configs_t* _search_instr;
    std::map<std::string, unsigned short> _typemap;
    std::vector< std::vector<std::string> > _combinations;
    CrystalSearch* _search;
    SEARCHTYPE Smode;
    BITTYPE Bmode;
    randomstore* _rdinfo;
    bool saveBest;
    double bestEnergy;
    
    
    double _density;
    double _lambda;
    unsigned int _iterations;
    unsigned int _popCount;
    unsigned int _nsweeps;
    
    unsigned int _maxbit;
    unsigned int _minbit;
    unsigned int _deltabit;

    
    virtual void _initialize();
    void keepBest(double);
    
    

    
    
};

#endif /* search_driver_hpp */

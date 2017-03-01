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

class SearchDriver{
public:
    
    enum SEARCHTYPE { ENANTIOPURE=0, RACEMIC=1, RACEMIZING=2};
    enum BITTYPE {RANDOMSWEEP=0, FIXEDSWEEP=1, NOSWEEP=2};
    
    SearchDriver(short, SEARCHTYPE,
                 BITTYPE,
                 double density,
                 double lambda,
                 int popCount,
                 int iterations);
    
    virtual void run();
    
    void setSearchType(SEARCHTYPE);
    void setBitType(BITTYPE);
    void setIterations(unsigned int);
    void setPopulationCount(unsigned int);
    //void setDensity(double);
    //void setTypeMap(std::map<std::string, unsigned short>&);
    //void setNumberOfMolecules(int);
    
    void overrideGenes(std::vector<double>&);
    
    
private:

protected:
    std::vector<molecule> _molecule_list;
    configs_t* _search_instr;
    std::map<std::string, unsigned short> _typemap;
    std::vector< std::vector<std::string> > _combinations;
    CrystalSearch* _search;
    SEARCHTYPE Smode;
    BITTYPE Bmode;
    
    double _density;
    double _lambda;
    unsigned int _iterations;
    unsigned int _popCount;
    
    unsigned int _maxbit;
    unsigned int _minbit;
    unsigned int _deltabit;
    
    virtual void _initialize();
    
    

    
    
};

#endif /* search_driver_hpp */

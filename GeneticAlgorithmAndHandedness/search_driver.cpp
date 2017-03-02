//
//  search_driver.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/23/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include "search_driver.hpp"
#include "molecules.h"

SearchDriver::SearchDriver(short nmol,
                           SEARCHTYPE st = ENANTIOPURE,
                           BITTYPE bt = RANDOMSWEEP,
                           double density= 0.24,
                           double lambda = 0.5,
                           int popCount= 10000,
                          int iterations= 10000):
_nmol(nmol),Smode(st),Bmode(bt),_density(density),
_lambda(lambda),_popCount(popCount),_iterations(iterations)
{
    _initialize();
}


void SearchDriver::_initialize(){

    _maxbit = 31;
    _minbit = 7;
    _deltabit = 2;

    _molecule_list.reserve(_nmol);
    
    if (Smode == ENANTIOPURE){
        for (auto& m : _molecule_list){
            m.second = "N";
            m.first = moleculestore.LTetramer;
        }
        std::vector< std::vector<std::string> > combns(1, std::vector<std::string>(2));
        combns[0][0] = "N"; combns[0][1] = "N";
        _combinations = combns;
        _typemap["N"] = 0;
    }
    else if (Smode == RACEMIC || Smode == RACEMIZING){
        int i=0;
        for (auto& m : _molecule_list){
            if (i%2 == 0){
                m.second = "N";
                m.first = moleculestore.LTetramer;
            }
            else{
                m.second = "O";
                m.first = moleculestore.DTetramer;
            }
        }
        std::vector< std::vector<std::string> > combns(1, std::vector<std::string>(2));
        combns[0][0] = "N"; combns[0][1] = "O";
        if (Smode == RACEMIZING){
            combns.push_back(combns[0]);
            combns[1][0] = "N"; combns[1][1] = "N";
        }
        _combinations = combns;
        _typemap["N"] = 0;
        _typemap["O"] = 1;
    }
    else{
        std::cerr << "Unkown switch mode.\n";
        exit(-1);
    }
    
    _search_instr = new configs_t(_nmol);
    _search_instr->setMolecules(_molecule_list);
    _search_instr->chiralitymap["N"] =  1.0;
     if (Smode == RACEMIZING || Smode == RACEMIC) _search_instr->chiralitymap["O"] = -1.0;
    _search_instr->populationCount=_popCount;
    _search_instr->iterations=_iterations;
    
    _search = new CrystalSearch(*_search_instr);
    _search->setDensity(_density);
    _search->setLambda(_lambda);
    _search->setTypeMap(_typemap);
    
    
}


void SearchDriver::setBitType(BITTYPE bt){
    Bmode = bt;
}

void SearchDriver::updateBittage(unsigned int bg){
    assert(bg <= 31);
    _search->updateBittage(bg);
}



/*class SearchDriver{
public:
    void run();
    
protected:

    
private:
    
};*/

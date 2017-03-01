//
//  search_driver.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/23/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include "search_driver.hpp"
#include "molecules.h"

void SearchDriver::_initialize(){
    _density = 0.24;
    _lambda = 0.50;
    _maxbit = 31;
    _minbit = 7;
    _deltabit = 2;
    _iterations = 10000;
    _popCount = 10000;
    
    _molecule_list.reserve(2);
    _molecule_list[0].second="N";
    _molecule_list[0].first = moleculestore.LTetramer;
    
    _molecule_list[1].second="O";
    _molecule_list[1].first = moleculestore.DTetramer;
    
    
    _combinations.reserve(2);
    _combinations[0].push_back("N");_combinations[0].push_back("N");
    _combinations[1].push_back("O");_combinations[1].push_back("O");
    
    _typemap["N"] = 0;
    _typemap["O"] = 1;
    
    _search_instr = new configs_t(2);
    _search_instr->setMolecules(_molecule_list);
    _search_instr->chiralitymap["N"] =  1.0;
    _search_instr->chiralitymap["O"] = -1.0;
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

void SearchDriver::setSearchType(SEARCHTYPE st){
    Smode = st;
}


/*class SearchDriver{
public:
    void run();
    
protected:

    
private:
    
};*/

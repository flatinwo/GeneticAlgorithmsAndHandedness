//
//  search_driver.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/23/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include "search_driver.hpp"
#include "molecules.h"
#include <limits>

static std::random_device rd;  //Will be used to obtain a seed for the random number engine
static std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

SearchDriver::SearchDriver(short nmol,
                           SEARCHTYPE st = ENANTIOPURE,
                           BITTYPE bt = RANDOMSWEEP,
                           double density= 0.24,
                           double lambda = 0.5,
                           int popCount= 10000,
                           int iterations= 100,
                           int bitmin=7,
                           int bitmax=31):
_nmol(nmol),Smode(st),Bmode(bt),_density(density),
_lambda(lambda),_popCount(popCount),_iterations(iterations),
_minbit(bitmin),_maxbit(bitmax)
{
    _search = nullptr;
    _search_instr = nullptr;
    _rdinfo = nullptr;
    _nsweeps=100;
    _initialize();
}

SearchDriver::~SearchDriver(){
    delete _search;
    delete _search_instr;
    if (_rdinfo != nullptr) delete _rdinfo;
}

void SearchDriver::setSeed(unsigned seed){
    assert(_rdinfo != nullptr);
    assert(0);
    
}

void SearchDriver::keepBest(double current_en){
   
    if (current_en < bestEnergy){
        std::ofstream logBest("logBest.dat",std::ios_base::app|std::ios_base::out);
        logBest << current_en << "\t" << bestEnergy << "\n";
        logBest.close();
        
        bestEnergy = current_en;
        _search->writeXYZ(2.0,"best_lattices.xyz");
        _search->printBestIndividual();
        
    }
    
}

void SearchDriver::run(){
    std::string str1("final_lattice_bit_");
    std::string str2(".xyz");
    std::string fname("");
    _search->writeXYZ(2.0,"initial_lattice.xyz");
    for (unsigned int i=0; i<_nsweeps;i++){
        
        if (Bmode == RANDOMSWEEP)
            _search->updateBittage(_rdinfo->generate());
        else if (Bmode == FIXEDSWEEP)
            _search->updateBittage(_minbit + i*_deltabit);
        else if (Bmode == NARROWINGRANDOMSWEEP)
            _search->updateBittage(_rdinfo->generate(_minbit+i*_deltabit,_maxbit));
        
        _search->run();
            
        fname = str1 + std::to_string(_search->getBittage())+ str2;
        _search->sortByFitness();
        _search->writeXYZ(2.0,fname.c_str());
        
        double current_energy = _search->getLatticeEnergy();
        if (saveBest) keepBest(current_energy);
        
        std::cout << "Current Best Lattice energy is: " << current_energy << "\n";
        if (saveBest) std::cout << "Overall Best Lattice energy is: " << bestEnergy << "\n";
    }
    
}

void SearchDriver::setIterationSweep(unsigned int sweeps ){
    _nsweeps = sweeps;
}

void SearchDriver::_initialize(){
    bestEnergy = std::numeric_limits<double>::max();
    saveBest = true;
    _deltabit = 1;
    if (Bmode == RANDOMSWEEP || Bmode == NARROWINGRANDOMSWEEP)
        _rdinfo = new randomstore(_minbit,_maxbit);

    
    _molecule_list.resize(_nmol,molecule());
    
    if (Smode == ENANTIOPURE){
        for (auto& m : _molecule_list){
            m.second = "N";
            m.first = moleculestore.LTetramer;
        }
        std::vector< std::vector<std::string> > combns(1, std::vector<std::string>(2));
        combns[0][0] = "N"; combns[0][1] = "N";
        if (_molecule_list.size() == 1) combns[0].pop_back();
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
            i++;
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
    _search_instr->setCombinations(_combinations);
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
    assert(bg <= 31 || bg >= 2);
    _search->updateBittage(bg);
    if (_rdinfo){
        if (bg >= _maxbit) _maxbit = bg;
        else _minbit = bg;
        _rdinfo->max = _maxbit;
        _rdinfo->min = _minbit;
    }
}



/*class SearchDriver{
public:
    void run();
    
protected:

    
private:
    
};*/

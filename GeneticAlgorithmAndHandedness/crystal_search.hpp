//
//  crystal_search.hpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 12/12/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef crystal_search_hpp
#define crystal_search_hpp

#include <stdio.h>
#include "operators.hpp"
#include "generation.hpp"
#include <string>
#include <fstream>
#include <array>
#include <algorithm>
#include <cassert>


#define MYCOMBMAX 32

typedef std::vector<gene> genes;

struct configs_t{
    std::vector<molecule> molecules;
    std::vector<gene> initGenes;
    
    short maxmolecules;
    short populationCount;
    short iterations;
    unsigned int bit;
    std::vector< std::vector<std::string> > combinations;
    unsigned long ncombinations=0;
    
    const bool o_flag;
    const bool s_flag;
    
    std::map< std::string, double> chiralitymap;
    
    configs_t(int n=1, bool or_flag=true, bool sw_flag=true, unsigned int xbit=10):bit(xbit),o_flag(or_flag),s_flag(sw_flag){
        molecules.resize(n,molecule());
        maxmolecules = 1024;
        iterations = 500;
        
        populationCount = 1000;
        
        //main lattice
        //idea make this a lattice object
        initGenes.push_back(gene(bit,0,1.0)); //x
        initGenes.push_back(gene(bit,0,1.0)); //y
        initGenes.push_back(gene(bit,0,PI/2.0)); //gamma
        initGenes.push_back(gene(bit,0,PI/2.0)); //phi
        initGenes.push_back(gene(bit,0,PI)); //psi
        
        for (unsigned int i=1; i<n; i++)
            for (unsigned int j=0;j<3; j++)
                initGenes.push_back(gene(bit,-1.0,1.0)); //cij -- for each center of mass displacement... change back to zero
        
        if (o_flag)
            for (unsigned int i=0; i<n; i++)
                for (unsigned int j=0;j<4;j++)
                    initGenes.push_back(gene(bit,-1.0,1.0)); //quaternions for rotation
        
        if (s_flag) initGenes.push_back(gene(2,0.0,1.0)); //probability of switching molecule identity
        
    }
    
    void setCombinations(const std::vector< std::vector<std::string> >& combns){
        assert(molecules.size() > 0);
        assert(s_flag);
        assert(combns.size() < MYCOMBMAX);
        
        unsigned int cbit = pow(2,combns.size()+1);
        for (auto const &n : combns) assert(n.size() < maxmolecules);
        
        initGenes.pop_back();
        initGenes.push_back(gene(cbit,0.0,1.0));
       
        //combinations.resize(combns.size());
        combinations = combns;
        ncombinations = combinations.size();
    }
    
    void setMolecules(const std::vector<molecule>& mols){
        molecules = mols;
    }
    
    void setIterations(unsigned int n){
        iterations = n;
    }
};

class CrystalSearch{
    friend class CrystalSearchTest_updateBittage_Test;
    friend class CrystalSearchTest_updateBittageValue_Test;
public:
    CrystalSearch();
    CrystalSearch(std::string);
    CrystalSearch(configs_t&);
    ~CrystalSearch();
    
    void writeXYZ(double=0.,const char* filename="lattice.xyz");
    void run();
    double getLatticeEnergy();
    
    void overrideGeneration(std::ifstream);
    void overrideGeneration();
    void overrideGeneration(std::vector<double>&);
    void sortByFitness();
    
    void setDensity(double);
    void setTypeMap(const std::map<std::string, unsigned short>&);
    void setLambda(double=0.);
    
    void updateBittage(unsigned int);
    unsigned getBittage();
    
    void resetIterationLog(){_itlog = 0;};
    void printBestIndividual();
    
    
protected:
    double _density;
    double _lambda;
    long _itlog;
    unsigned int _gencount;
    generation* _generation;
    const configs_t*  _myconfig;
    individual* _ind;
    genes* _genes;
    lattice _lattice;
    
    double3 _a1,_a2,_a3; // lattice
    
    std::vector<double3> _dcm; // displacement of center of mass
    std::vector<double4> _quat; // quaternion vector
    std::vector<double4> _axangle; //axis angle representation
    
    
    std::map<std::string,unsigned short> _typemap;
    
    std::vector<individual> _individuals;
    std::vector<molecule> _mymolecules;
    std::map< std::string, double> _chiralitymap;
    
    ofstream* fitInfo;
    ofstream* solPInfo;
    ofstream* solVInfo;
    
    void _initialize();
    void _openFiles();
    void _closeFiles();
    void _printIndividual(unsigned int);
    void _writeXYZ(double=0., const char* filename="lattice.xyz");
    
    double _computeFitness(individual&, bool returnEnergy=false); // is static function inherited
    void _computeGenerationFitness();
    void _setUp(individual*,bool=false); //to set up values for fitness and writing
    
    enum LOCALINFO {ROTATE=0, TRANSLATE=1, SDIM=3, ODIM=4};
    
    std::array<double,LOCALINFO::SDIM> _cms;
    std::array<double,LOCALINFO::ODIM> _qts;
    
    
    
};

#endif /* crystal_search_hpp */

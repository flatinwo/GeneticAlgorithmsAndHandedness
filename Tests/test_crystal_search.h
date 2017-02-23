//
//  test_crystal_search.h
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/23/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#ifndef test_crystal_search_h
#define test_crystal_search_h


#include "crystal_search.hpp"
#include "gtest/gtest.h"
#include "GenAlgoAndHands.h"

struct CrystalSearchTest : public testing::Test{
    CrystalSearch* search;
    int iterations=20;
    int popCount=20;
    double density=0.2;
    CrystalSearchTest(){
        const int nmol = 2;
        
        std::vector<molecule> twomolecules(nmol);
        twomolecules[0].second = "N";
        twomolecules[1].second = "O";
        
        std::vector<double3> refTetramer(4);
        refTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
        refTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
        refTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
        refTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
        
        refTetramer[0]=double3(3.41961,5.16472,10.6341);
        refTetramer[1]=double3(3.22152,6.20496,10.6203);
        refTetramer[2]=double3(4.26042,6.40321,10.5689);
        refTetramer[3]=double3(4.18845,6.51404,9.51877);
        
        
        twomolecules[0].first = twomolecules[1].first = refTetramer;
        twomolecules[1].first[3] = double3( 0.52915050,  0.26457525, 1.32287625);
        
        twomolecules[1].first[0] = double3( 3.44621, 4.62218, 9.61783);
        twomolecules[1].first[1] = double3( 3.56946, 5.66827, 9.71798);
        twomolecules[1].first[2] = double3( 4.60823, 5.53141, 9.86875);
        twomolecules[1].first[3] = double3( 4.44067, 5.39446, 10.9049);
        
        std::vector< std::vector<std::string> > combns(1, std::vector<std::string>(2));
        combns[0][0] = "N"; combns[0][1] = "O";
        
        std::map<std::string,unsigned short> typemap;
        typemap["N"] = 0; typemap["O"] = 1;
        

        
        
        configs_t search_instr(nmol); //search instructions (2 molecules per search)
        search_instr.setMolecules(twomolecules);
        search_instr.chiralitymap["N"] =  1.0;
        search_instr.chiralitymap["O"] =  -1.0;
        search_instr.setCombinations(combns);
        search_instr.iterations = iterations;
        search_instr.populationCount=popCount;
        
        search = new CrystalSearch(search_instr);
        search->setDensity(density);
        search->setTypeMap(typemap);
    }
    
    virtual ~CrystalSearchTest(){
        delete search;
    }
};

TEST_F(CrystalSearchTest,updateBittage){
    for (int i=0 ; i<popCount; i++){
        genes* gn = &(search->_generation->individuals[i].genes);
        for (std::vector<gene>::iterator it = gn->begin();
             it != gn->end()-1;  it++) EXPECT_EQ(10,it->code.size());
        EXPECT_EQ(4,gn->back().code.size());
    }
    
    //now change bittage
    search->updateBittage(14);
    
    for (int i=0 ; i<popCount; i++){
        genes* gn = &(search->_generation->individuals[i].genes);
        for (std::vector<gene>::iterator it = gn->begin();
             it != gn->end()-1;  it++) EXPECT_EQ(14,it->code.size());
        EXPECT_EQ(4,gn->back().code.size());
    }
    
}

TEST_F(CrystalSearchTest,updateBittageValue){
    
}

#endif /* test_crystal_search_h */

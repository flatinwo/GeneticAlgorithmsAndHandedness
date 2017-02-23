//
//  test_gene.h
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/6/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

//http://www.yolinux.com/TUTORIALS/Cpp-GoogleTest.html
//https://blog.jetbrains.com/rscpp/unit-testing-google-test/

#ifndef test_gene_h
#define test_gene_h

#include "gene.hpp"
#include "gtest/gtest.h"
#include "GenAlgoAndHands.h"

struct geneTest : public testing::Test{
    gene* Gene;
    geneTest(){
        Gene = new gene(10,0.,9.);
    }
    
    virtual ~geneTest(){
        delete Gene;
    }
};

TEST_F(geneTest,decode){
    const double three=0.5;
    Gene->setValue(three);
    EXPECT_EQ(10,Gene->code.size());
    EXPECT_NEAR(0.5,Gene->getValue(),0.1);
    
}

TEST_F(geneTest,falsestrings){
    const double three=0.5;
    Gene->setValue(three);
    EXPECT_EQ(10,Gene->code.size());
    EXPECT_STRNE(NULL,Gene->makeString().c_str());
    
}

TEST_F(geneTest,BoundaryValues){
    const double four=8.;
    Gene->setValue(four);
    EXPECT_EQ(0.,Gene->getMinValue());
    EXPECT_EQ(9.,Gene->getMaxValue());
    EXPECT_NEAR(9./(pow(2.,10)),Gene->getResolution(),0.001);
    
}

TEST_F(geneTest,UpdateBittage){
    const double four=8.;
    Gene->setValue(four);
    EXPECT_NEAR(9./(pow(2.,10)),Gene->getResolution(),0.001);
    EXPECT_NEAR(8.,Gene->getValue(), Gene->getResolution());
    
    double value = Gene->getValue();
    unsigned int newbit=31;
    Gene->updateBittage(newbit);
    EXPECT_EQ(newbit,Gene->code.size());
    EXPECT_NEAR(9./(pow(2.,newbit)),Gene->getResolution(),0.001);
    EXPECT_NEAR(value, Gene->getValue(), Gene->getResolution());
    
    
}


#endif /* test_gene_h */

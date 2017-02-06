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
        Gene = new gene(10.,0.,9.);
    }
    
    virtual ~geneTest(){
        delete Gene;
    }
};

TEST_F(geneTest,decode){
    const double three=3.;
    Gene->setValue(three);
    EXPECT_EQ(10,Gene->code.size());
    EXPECT_NEAR(3.,Gene->getValue(),0.1);
    
}

TEST_F(geneTest,falsestrings){
    const double three=3.;
    Gene->setValue(three);
    EXPECT_EQ(10,Gene->code.size());
    EXPECT_STRNE(NULL,Gene->makeString().c_str());
    
}



#endif /* test_gene_h */

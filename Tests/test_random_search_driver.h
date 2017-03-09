//
//  test_random_search_driver.h
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 3/6/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#ifndef test_random_search_driver_h
#define test_random_search_driver_h

#include "search_driver.hpp"
#include "gtest/gtest.h"
#include "GenAlgoAndHands.h"
#include <list>

struct TestRandomEngineSearchDriver: public testing::Test {
    randomstore* rst;
    double exp_u_mean;
    double exp_u_std;
    TestRandomEngineSearchDriver(unsigned a=0, unsigned b=10){
        rst = new randomstore(a,b);
        exp_u_mean = 0.5*(a+b);
        exp_u_std = sqrt((b-a)*(b-a)/12.);
    }
    
    ~TestRandomEngineSearchDriver(){
        delete rst;
    }
    
};

TEST_F(TestRandomEngineSearchDriver, checkMean){
    std::list<unsigned> samplesizes{500,1000,5000};
    for (std::list<unsigned>::iterator it=samplesizes.begin();
         it != samplesizes.end(); it++){
        float sum=0; int i=0;
        while (i < *it){
            sum += rst->generate();
            i++;
        }
        double mean = sum/double(*it);
        double conf_interval = 1.96*exp_u_std/sqrt(*it);
        EXPECT_NEAR(mean,exp_u_mean, conf_interval);
        
    }
    
}


#endif /* test_random_search_driver_h */

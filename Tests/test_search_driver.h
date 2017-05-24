//
//  test_search_driver.h
//  GeneticAlgorithmAndHandedness
//
//  Created by DrFlo on 5/24/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#ifndef test_search_driver_h
#define test_search_driver_h

#include "search_driver.hpp"
#include "gtest/gtest.h"
#include "GenAlgoAndHands.h"

struct SearchDriverTest: public testing::Test{
    SearchDriver* driver;
    SearchDriverTest(){};
    virtual ~SearchDriverTest(){};
    
    void SetUp(short ind, double density, double alpha){
        switch (ind) {
            case 0:
            {
                const short n = 1;
                driver = new  SearchDriver(n,
                                           SearchDriver::SEARCHTYPE::ENANTIOPURE,
                                           SearchDriver::BITTYPE::RANDOMSWEEP,
                                           density,
                                           alpha,
                                           10,
                                           2,
                                           31,
                                           31);
                break;
            }

            case 1:
            {
                const short n = 2;
                driver = new  SearchDriver(n,
                                           SearchDriver::SEARCHTYPE::RACEMIZING,
                                           SearchDriver::BITTYPE::RANDOMSWEEP,
                                           density,
                                           alpha,
                                           10,
                                           2,
                                           31,
                                           31);
                break;
            }

            case 2:
            {
                const short n = 2;
                driver = new  SearchDriver(n,
                                           SearchDriver::SEARCHTYPE::RACEMIC,
                                           SearchDriver::BITTYPE::RANDOMSWEEP,
                                           density,
                                           alpha,
                                           10,
                                           2,
                                           31,
                                           31);
                driver->updateBittage(17);
                break;
            }
            case 4:
            {
                const short n = 2;
                driver = new  SearchDriver(n,
                                           SearchDriver::SEARCHTYPE::ENANTIOPURE,
                                           SearchDriver::BITTYPE::RANDOMSWEEP,
                                           density,
                                           alpha,
                                           10,
                                           2,
                                           31,
                                           31);
                break;
            }
                
            default:
                std::cerr <<"Unknown choice for switch option\n";
                break;
        }
    }
    
    
    virtual void TearDown(){
        delete driver;
    }
};

TEST_F(SearchDriverTest,Energy1){
    
    SetUp(2,0.24,-0.50);
    std::vector<double> v{0.789666,	0.624805,	0.508338,	0.665289,	0.6255,	-0.0182806,	0.124851,	0.499619,	-0.0159612,	0.851497,	0.0168157,	-0.832514,	-0.562731,	-0.575915,	-0.600208,	0.588855,	0.733333};
    driver->overrideGenes(v);
   EXPECT_NE(-30.3527, driver->getEnergy()) ;
}


#endif /* test_search_driver_h */

//
//  main.cpp
//  Tests
//
//  Created by Folarin Latinwo on 2/6/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include <iostream>
#include "gtest/gtest.h"
#include "test_crystal_search.h"
#include "test_lattice.h"
#include "test_gene.h"
#include "gene.hpp"
#include "test_random_search_driver.h"
#include "test_search_driver.h"

int main(int argc, char * argv[]) {

    
    std::cout << "Hello, World!\n";

    
    // insert code here...
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

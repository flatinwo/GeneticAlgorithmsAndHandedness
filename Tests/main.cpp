//
//  main.cpp
//  Tests
//
//  Created by Folarin Latinwo on 2/6/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include <iostream>
#include "gtest/gtest.h"
#include "test_gene.h"
#include "gene.hpp"

int main(int argc, char * argv[]) {

    
    std::cout << "Hello, World!\n";

    
    // insert code here...
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

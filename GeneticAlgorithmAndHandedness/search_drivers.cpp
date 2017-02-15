//
//  search_drivers.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 2/6/17.
//  Copyright © 2017 Folarin Latinwo. All rights reserved.
//

#include "search_drivers.hpp"
#include "crystal_search.hpp"


void TestNmolecule1(double density=0.24, double alpha=0.5,
                    int iterations=2000, int popCount=5000){
    const int nmol = 1;
    
    std::vector<molecule> onemolecules(nmol);
    onemolecules[0].second = "N";
    
    std::vector<double3> refTetramer(4);
    refTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
    refTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
    refTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
    refTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
    
    
    onemolecules[0].first = refTetramer;
    
    std::vector< std::vector<std::string> > combns(1, std::vector<std::string>(1)); //can implement case for switching number of basis pts.
    combns[0][0] = "N";
    
    std::map<std::string,unsigned short> typemap;
    typemap["N"] = 0;
    
    
    configs_t search_instr(nmol); //search instructions (2 molecules per search)
    search_instr.setMolecules(onemolecules);
    search_instr.chiralitymap["N"] =  1.0;
    search_instr.setCombinations(combns);
    search_instr.populationCount=popCount;
    search_instr.iterations = iterations;
    
    CrystalSearch onemolsearch(search_instr);
    onemolsearch.setDensity(density);
    
    /*std::vector<double> overrideGenes(10,0.);
     overrideGenes[ 0]=1.0;overrideGenes[ 1]=1.0;overrideGenes[ 2]=PI/2;overrideGenes[ 3]=0.0;overrideGenes[ 4]=PI/2;
     overrideGenes[ 5]=1.0;overrideGenes[ 6]=0.0;overrideGenes[ 7]=0.0;overrideGenes[ 8]=0.0;
     overrideGenes[9]=0.1;
     
     
     overrideGenes[ 0]=0.747986;overrideGenes[ 1]=0.982421;overrideGenes[ 2]=1.21303;overrideGenes[ 3]=0.197094;overrideGenes[ 4]=1.565;
     overrideGenes[ 5]=0.843749;overrideGenes[ 6]=-0.187775;overrideGenes[ 7]=-0.0717774;overrideGenes[ 8]=0.358588;
     overrideGenes[9]=0.0531288;
     
     onemolsearch.overrideGeneration(overrideGenes);*/
    
    
    onemolsearch.setLambda(alpha);
    onemolsearch.setTypeMap(typemap);
    
    onemolsearch.writeXYZ(2.0);
    std::cout << "Lattice energy is: " << onemolsearch.getLatticeEnergy() << "\n";
    
    onemolsearch.run();
    
    std::cout << "Lattice energy is: " << onemolsearch.getLatticeEnergy() << "\n";
    onemolsearch.writeXYZ(2.0);
    
}





void TestNmolecule2(double density=0.24, double alpha=0.5,
                    int iterations=2000, int popCount=5000){
    const int nmol = 2;
    
    std::vector<molecule> twomolecules(nmol);
    twomolecules[0].second = "N";
    twomolecules[1].second = "O";
    
    std::vector<double3> refTetramer(4);
    refTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
    refTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
    refTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
    refTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
    
    
    twomolecules[0].first = twomolecules[1].first = refTetramer;
    twomolecules[1].first[3] = double3( 0.52915050,  0.26457525, 1.32287625);
    
    std::vector< std::vector<std::string> > combns(2, std::vector<std::string>(2)); //can implement case for switching number of basis pts.
    combns[0][0] = "N"; combns[0][1] = "N";
    combns[1][0] = "N"; combns[1][1] = "O";
    
    std::map<std::string,unsigned short> typemap;
    typemap["N"] = 0; typemap["O"] = 1;
    
    
    configs_t search_instr(nmol); //search instructions (2 molecules per search)
    search_instr.setMolecules(twomolecules);
    search_instr.chiralitymap["N"] =  1.0;
    search_instr.chiralitymap["O"] =  -1.0;
    search_instr.setCombinations(combns);
    search_instr.populationCount=popCount;
    search_instr.iterations = iterations;
    
    CrystalSearch twomolsearch(search_instr);
    twomolsearch.setDensity(density);
    
    /*std::vector<double> overrideGenes(17,0.);
     overrideGenes[ 0]=1.0;overrideGenes[ 1]=1.0;overrideGenes[ 2]=PI/2;overrideGenes[ 3]=0.0;overrideGenes[ 4]=PI/2;
     overrideGenes[ 5]=0.22;overrideGenes[ 6]=0.3;overrideGenes[ 7]=0.3 ;
     overrideGenes[ 8]=1.0;overrideGenes[ 9]=0.0;overrideGenes[10]=0.0;overrideGenes[11]=0.0;
     overrideGenes[12]=1.0;overrideGenes[13]=0.0;overrideGenes[14]=0.0;overrideGenes[15]=0.0;
     overrideGenes[16]=0.6;*/
    
    //twomolsearch.overrideGeneration(overrideGenes);
    
    
    twomolsearch.setLambda(alpha);
    twomolsearch.setTypeMap(typemap);
    twomolsearch.writeXYZ(2.0);
    std::cout << "Lattice energy is: " << twomolsearch.getLatticeEnergy() << "\n";
    
    
    twomolsearch.run();
    
    twomolsearch.writeXYZ(2.0);
    std::cout << "Lattice energy is: " << twomolsearch.getLatticeEnergy() << "\n";
    
    
}

void TestNmolecule2E(double density=0.24, double alpha=0.5,
                     int iterations=2000, int popCount=5000){
    const int nmol = 2;
    
    std::vector<molecule> twomolecules(nmol);
    twomolecules[0].second = "N";
    twomolecules[1].second = "N";
    
    std::vector<double3> refTetramer(4);
    refTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
    refTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
    refTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
    refTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
    
    
    
    twomolecules[0].first = twomolecules[1].first = refTetramer;

    
    std::vector< std::vector<std::string> > combns(1, std::vector<std::string>(2)); //can implement case for switching number of basis pts.
    combns[0][0] = "N"; combns[0][1] = "N";
    //combns[1][0] = "N"; combns[1][1] = "O";
    
    std::map<std::string,unsigned short> typemap;
    typemap["N"] = 0; typemap["O"] = 1;
    
    
    configs_t search_instr(nmol); //search instructions (2 molecules per search)
    search_instr.setMolecules(twomolecules);
    search_instr.chiralitymap["N"] =  1.0;
    //search_instr.chiralitymap["O"] =  -1.0;
    search_instr.setCombinations(combns);
    search_instr.iterations = iterations;
    search_instr.populationCount=popCount;
    
    CrystalSearch twomolsearch(search_instr);
    twomolsearch.setDensity(density);
    twomolsearch.setTypeMap(typemap);
    
    //std::vector<double> myGenes({0.500023,	0.750408,	1.37391,	0.49233,	0.585289,	0.822215,	0.864837,	0.629519,	-0.00186597,	0.962939,	-0.764922, -0.625887,	0.50115,	0.0719044,	-0.998391, -0.20599,	0.560162});
    
    //std::vector<double> myGenes({1.,	1.,	PI/2,	0.,	PI/2,	0.2436425,	-0.7676525,	-0.3081525,	1.,	.0,	0., 0.,	1.,	.0,	0., 0.,	0.1});
    
    /*std::vector<double> myGenes({ 0.497024,	0.928539,	0.676501,	0.556789,	0.882762,	-0.374863,	-0.375708,	0.685241,	0.796704,	0.148485,	0.755849,	-0.218731,	0.561221,	-0.438728,	-0.647232,	-0.584263,	0.753198});*/
    
    
    /*std::vector<double> overrideGenes(17,0.);
     overrideGenes[ 0]=1.0;overrideGenes[ 1]=1.0;overrideGenes[ 2]=PI/2;overrideGenes[ 3]=0.0;overrideGenes[ 4]=PI/2;
     overrideGenes[ 5]=0.22;overrideGenes[ 6]=0.3;overrideGenes[ 7]=0.3 ;
     overrideGenes[ 8]=1.0;overrideGenes[ 9]=0.0;overrideGenes[10]=0.0;overrideGenes[11]=0.0;
     overrideGenes[12]=1.0;overrideGenes[13]=0.0;overrideGenes[14]=0.0;overrideGenes[15]=0.0;
     overrideGenes[16]=0.6;*/
    
    
    //twomolsearch.overrideGeneration(myGenes);
    twomolsearch.setLambda(alpha);
    
    std::cout << "Lattice energy is: " << twomolsearch.getLatticeEnergy() << "\n";
    twomolsearch.writeXYZ(2.0);
    
    
    twomolsearch.run();
    
    twomolsearch.writeXYZ(2.0);
    std::cout << "Lattice energy is: " << twomolsearch.getLatticeEnergy() << "\n";
    
    
}



void TestNmolecule2R(double density=0.24, double alpha=0.5,
                     int iterations=2000, int popCount=5000){
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
    
    std::vector< std::vector<std::string> > combns(1, std::vector<std::string>(2)); //can implement case for switching number of basis pts.
    combns[0][0] = "N"; combns[0][1] = "O";
    //combns[1][0] = "N"; combns[1][1] = "O";
    
    std::map<std::string,unsigned short> typemap;
    typemap["N"] = 0; typemap["O"] = 1;
    
    
    configs_t search_instr(nmol); //search instructions (2 molecules per search)
    search_instr.setMolecules(twomolecules);
    search_instr.chiralitymap["N"] =  1.0;
    search_instr.chiralitymap["O"] =  -1.0;
    search_instr.setCombinations(combns);
    search_instr.iterations = iterations;
    search_instr.populationCount=popCount;
    
    CrystalSearch twomolsearch(search_instr);
    twomolsearch.setDensity(density);
    twomolsearch.setTypeMap(typemap);
    
    //std::vector<double> myGenes({0.500023,	0.750408,	1.37391,	0.49233,	0.585289,	0.822215,	0.864837,	0.629519,	-0.00186597,	0.962939,	-0.764922, -0.625887,	0.50115,	0.0719044,	-0.998391, -0.20599,	0.560162});
    
    //std::vector<double> myGenes({1.,	1.,	PI/2,	0.,	PI/2,	0.2436425,	-0.7676525,	-0.3081525,	1.,	.0,	0., 0.,	1.,	.0,	0., 0.,	0.1});
    
    std::vector<double> myGenes({ 0.497024,	0.928539,	0.676501,	0.556789,	0.882762,	-0.374863,	-0.375708,	0.685241,	0.796704,	0.148485,	0.755849,	-0.218731,	0.561221,	-0.438728,	-0.647232,	-0.584263,	0.753198});
    
    
    /*std::vector<double> overrideGenes(17,0.);
     overrideGenes[ 0]=1.0;overrideGenes[ 1]=1.0;overrideGenes[ 2]=PI/2;overrideGenes[ 3]=0.0;overrideGenes[ 4]=PI/2;
     overrideGenes[ 5]=0.22;overrideGenes[ 6]=0.3;overrideGenes[ 7]=0.3 ;
     overrideGenes[ 8]=1.0;overrideGenes[ 9]=0.0;overrideGenes[10]=0.0;overrideGenes[11]=0.0;
     overrideGenes[12]=1.0;overrideGenes[13]=0.0;overrideGenes[14]=0.0;overrideGenes[15]=0.0;
     overrideGenes[16]=0.6;*/
    
    
    twomolsearch.overrideGeneration(myGenes);
    twomolsearch.setLambda(alpha);
    
    std::cout << "Lattice energy is: " << twomolsearch.getLatticeEnergy() << "\n";
    twomolsearch.writeXYZ(2.0);
    
    
    twomolsearch.run();
    
    twomolsearch.writeXYZ(2.0);
    std::cout << "Lattice energy is: " << twomolsearch.getLatticeEnergy() << "\n";
    
    
}

void TestNmolecule4R(double density=0.24, double alpha=0.5,
                     int iterations=2000, int popCount=5000){
    const int nmol = 4;
    
    std::vector<molecule> fourmolecules(nmol);
    fourmolecules[0].second = "N1";
    fourmolecules[1].second = "O1";
    fourmolecules[2].second = "N2";
    fourmolecules[3].second = "O2";
    
    std::vector<double3> refTetramer(4);
    refTetramer[0] = double3(-0.52915050, -0.79372575,  0.26457525);
    refTetramer[1] = double3(-0.52915050,  0.26457525,  0.26457525);
    refTetramer[2] = double3( 0.52915050,  0.26457525,  0.26457525);
    refTetramer[3] = double3( 0.52915050,  0.26457525, -0.79372575);
    
    
    fourmolecules[0].first = fourmolecules[1].first = fourmolecules[2].first = fourmolecules[3].first = refTetramer;
    fourmolecules[1].first[3] = fourmolecules[3].first[3] = double3( 0.52915050,  0.26457525, 1.32287625);
    
    std::vector< std::vector<std::string> > combns(1, std::vector<std::string>(4)); //can implement case for switching number of basis pts.
    combns[0] = {"N1","O1","N2","O2"};
    //combns[1][0] = "N"; combns[1][1] = "O";
    
    std::map<std::string,unsigned short> typemap;
    typemap["N1"] = 0; typemap["O1"] = 1;
    typemap["N2"] = 2; typemap["O2"] = 3;
    
    
    configs_t search_instr(nmol); //search instructions (2 molecules per search)
    search_instr.setMolecules(fourmolecules);
    search_instr.chiralitymap["N1"] = search_instr.chiralitymap["N2"] = 1.0;
    search_instr.chiralitymap["O1"] = search_instr.chiralitymap["O2"] = -1.0;
    search_instr.setCombinations(combns);
    search_instr.iterations = iterations;
    search_instr.populationCount=popCount;
    
    CrystalSearch fourmolsearch(search_instr);
    fourmolsearch.setDensity(density);
    fourmolsearch.setTypeMap(typemap);
    
    //std::vector<double> myGenes({0.608201,	0.994877,	0.121688,	0.161599,	2.88972,	0.885562,	0.402058,	0.909395,	0.713875,	0.439257,	0.942913,	0.678971,	0.817784,	0.852869,	0.799386,	0.876072,	-0.769296,	-0.412556,	0.112058,	0.137701,	-0.744683,	0.969938,	0.70773,	0.823122,	0.887617,	0.262039,	0.849356,	0.724021,	-0.473222,	0.813087,	0.321229});
    
    /*std::vector<double> overrideGenes(17,0.);
     overrideGenes[ 0]=1.0;overrideGenes[ 1]=1.0;overrideGenes[ 2]=PI/2;overrideGenes[ 3]=0.0;overrideGenes[ 4]=PI/2;
     overrideGenes[ 5]=0.22;overrideGenes[ 6]=0.3;overrideGenes[ 7]=0.3 ;
     overrideGenes[ 8]=1.0;overrideGenes[ 9]=0.0;overrideGenes[10]=0.0;overrideGenes[11]=0.0;
     overrideGenes[12]=1.0;overrideGenes[13]=0.0;overrideGenes[14]=0.0;overrideGenes[15]=0.0;
     overrideGenes[16]=0.6;*/
    
    
    //fourmolsearch.overrideGeneration(myGenes);
    fourmolsearch.setLambda(alpha);
    
    std::cout << "Lattice energy is: " << fourmolsearch.getLatticeEnergy() << "\n";
    fourmolsearch.writeXYZ(2.0);
    
    
    fourmolsearch.run();
    
    fourmolsearch.writeXYZ(2.0);
    std::cout << "Lattice energy is: " << fourmolsearch.getLatticeEnergy() << "\n";
    
    
}

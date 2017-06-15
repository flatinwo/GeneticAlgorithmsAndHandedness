//
//  crystal_search.cpp
//  GeneticAlgorithmAndHandedness
//
//  Created by Folarin Latinwo on 12/12/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "crystal_search.hpp"
#include <algorithm>
#include "operators.hpp"

#define ESPILON 0.0001
#define ONEMINUSEPSILON 0.9999

CrystalSearch::CrystalSearch(){
    _myconfig = nullptr;
    _initialize();
}

CrystalSearch::CrystalSearch(configs_t& cfig_t){
    _myconfig = &cfig_t;
    _initialize();
}

CrystalSearch::~CrystalSearch(){
    delete _generation;
    _closeFiles();
}

void CrystalSearch::setDensity(double dens){
    _density = dens;
}

void CrystalSearch::setLambda(double lam){
    _lambda = lam;
}

void CrystalSearch::overrideGeneration(){
    string input_line;
    
    std::vector<gene> iGenes(_myconfig->initGenes);
    unsigned count = 0;
    while (std::cin){
        getline(std::cin, input_line);
        iGenes[count++].setValue(std::stod(input_line)); //need some checks for this
        if (count > iGenes.size()) break;
    }
    
    for (unsigned int j=0; j<_generation->individuals.size();j++)
        _generation->individuals[j] = iGenes;
}

void CrystalSearch::writeXYZ(double myrcut,const char* filename){
    _writeXYZ(myrcut,filename);
}

void CrystalSearch::overrideGeneration(std::vector<double>& iGenes){
    assert(iGenes.size()  == _myconfig->initGenes.size());
    for (unsigned int j=0 ; j <  _generation->individuals.size() ; j++)
        for (unsigned int k=0; k<_generation->individuals[j].genes.size(); k++)
            if (!(_generation->individuals[j].genes[k].setValue(iGenes[k]))) exit(-1);
    
}

//note changing myconfig genes so as to have memory of initial bittage.
void CrystalSearch::updateBittage(unsigned int length){
    for (unsigned int i=0; i < _generation->individuals.size(); i++)
        for_each(_generation->individuals[i].genes.begin(),
                 _generation->individuals[i].genes.end()-1,
                 [length](gene& m){ m.updateBittage(length);}); //exclude orientation.
    
    
}

unsigned CrystalSearch::getBittage(){
    return (unsigned)_generation->individuals[0].genes[0].code.size();
}
void CrystalSearch::setTypeMap(const std::map<std::string, unsigned short>& tmap){
    _typemap = tmap;
}

void CrystalSearch::_initialize(){
    if (_myconfig == nullptr) _myconfig = new configs_t(2,true,true);
    _individuals.resize(_myconfig->populationCount);
    _mymolecules = _myconfig->molecules;
    
    std::fill(_individuals.begin(),_individuals.end(),_myconfig->initGenes);
    
    _generation = new generation(_individuals);
    _generation->randomize();
    _generation->elitism = true;
    
    _dcm.resize(_myconfig->molecules.size(),double3(0.,0.,0.));
    if (_myconfig->o_flag) {
        _quat.resize(_myconfig->molecules.size());
        _axangle.resize(_myconfig->molecules.size());
    }
    
    _gencount = 0;
    _openFiles();
    _chiralitymap = _myconfig->chiralitymap;
    _itlog = 0;
}

double CrystalSearch::getLatticeEnergy(){
    return _computeFitness(_generation->individuals[0],true);
}

double CrystalSearch::getClusterEnergy(){
    double en = 0.;
    _setUp(&(_generation->individuals[0]),true);
    
    for (unsigned int l=0;l<_mymolecules.size(); l++){
        double zeta1 = _chiralitymap[_mymolecules[l].second];
        for (unsigned int m=0; m<_mymolecules.size(); m++){
            if (l==m) continue;
            double zeta2 = _chiralitymap[_mymolecules[m].second];
            
            for (unsigned int n=0; n < _mymolecules[l].first.size(); n++){
                double3 ref_part = _mymolecules[l].first[n];
                for (unsigned int o=0; o<_mymolecules[m].first.size(); o++){
                    double3 part = _mymolecules[m].first[o];
                    en += (1. + _lambda*zeta1*zeta2)*vLJ(1./((part-ref_part).get2Norm()));
                }
            }
            
        }
    }
    return en/(double) _mymolecules.size();
}

double CrystalSearch::_computeFitness(individual& ind, bool returnEnergy){
    _setUp(&ind,true);
    double a1l, a2l, a3l;
    const double myrCut = 4.0;
    a1l = _a1.getLength(); a2l = _a2.getLength(); a3l = _a3.getLength();
    
    const int iMax = ceil(myrCut/a1l) + 1;
    const int jMax = ceil(myrCut/a2l) + 1;
    const int kMax = ceil(myrCut/a3l) + 1;
    
    double energy=0.;
    for (int i=-iMax; i<= iMax; i++){
        for (int j=-jMax; j<= jMax; j++ ){
            for (int k=-kMax; k<=kMax; k++){
                const double3 site_ijk = double3(i*_a1.x + j*_a2.x + k*_a3.x,
                                                 i*_a1.y + j*_a2.y + k*_a3.y,
                                                 i*_a1.z + j*_a2.z + k*_a3.z);
                
                for (unsigned int l=0;l<_mymolecules.size(); l++){
                    double zeta1 = _chiralitymap[_mymolecules[l].second];
                    
                    for (unsigned int m=0; m<_mymolecules.size(); m++){
                        if ((i==j) && (j==k) && (k==0) && l == m) continue;
                        
                        double zeta2 = _chiralitymap[_mymolecules[m].second];
                        //std::cout << "zeta1 and zeta2 are " << zeta1 << "\t" << zeta2 << "\n\n";
                        
                        for (unsigned int n=0; n < _mymolecules[l].first.size(); n++){
                            double3 ref_part = _mymolecules[l].first[n];
                            for (unsigned int o=0; o<_mymolecules[m].first.size(); o++){
                                double3 part = _mymolecules[m].first[o] + site_ijk;
                                energy += (1. + _lambda*zeta1*zeta2)*vLJ(1./((part-ref_part).get2Norm()));
                            }
                        }
                        
                    }
                    
                }
            }
        }
    }
    
    energy *= 0.5; //to avoid double counting
    energy /= _mymolecules.size();
    
    ind.fitness = exp(-1.*energy);
    
    if (returnEnergy) return energy;
    else return exp(-1.*energy);
}


void CrystalSearch::_closeFiles(){
    fitInfo->close();
    solPInfo->close();
    solVInfo->close();
    
    delete fitInfo;
    delete solPInfo;
    delete solVInfo;
}

void CrystalSearch::_openFiles(){
    fitInfo = new ofstream("fitness.dat");
    solPInfo = new ofstream("solutionParameter.dat");
    solVInfo = new ofstream("solutionVector.dat");
    
    //Information for solution parameters
    *solPInfo << "#g\t";
    for (unsigned int i=0; i < _myconfig->initGenes.size(); i++) *solPInfo << i << "\t";
    *solPInfo << "\n";
    
    //Information for solution vector
    *solVInfo << "#g\t";
    unsigned int count=9; //3*3
    count += 3*(_myconfig->molecules.size()); //3 center of mass coordinate
    if (_myconfig->o_flag) count += 4*(_myconfig->molecules.size()); //quaternions
    if (_myconfig->s_flag) count += 1;
    for (unsigned int i=0; i < count; i++) *solVInfo << i << "\t";
    *solVInfo << "\n";
    
    *fitInfo << "#g\tfitness\n";
    
}

void CrystalSearch::_computeGenerationFitness(){
    double fit=0.;
    for (auto& ind : _generation->individuals) fit = _computeFitness(ind);
    
    _generation->renormalizeFitness();
}

void CrystalSearch::run(){
    for (unsigned int t=0; t<_myconfig->iterations; t++){
        //set up generation
        //_generation->calcFitness(_computeFitness);
        _computeGenerationFitness();
        _generation->sortByFitness();
        
        //output generation
        *fitInfo << _itlog << "\t" << _generation->individuals[0].fitness <<
        "\t" << _generation->individuals[0].genes[0].code.size() << "\n";
        
        std::cout << _itlog << "\t" << _generation->individuals[0].fitness<<
        "\t" << _generation->individuals[0].genes[0].code.size()<< "\n";
        _setUp(&(_generation->individuals[0]),true);
        
        std::cout << "Fate: ";
        for (auto& a : _mymolecules) std::cout << a.second << "\t";
        std::cout << "\n";
        _printIndividual(0);
        
        //create next generation
        _generation->mateFitnessProportional();
        _generation->mutate(1.0, 0.05); //set parameters
        
        _gencount++;
        _itlog++;
    }
    
    _computeGenerationFitness();
    _generation->sortByFitness();
    
    std::cout << "Best fitness (after hill-search) is " << _computeFitness(_generation->individuals[0]) << "\n";
   // _writeXYZ();
}

void CrystalSearch::sortByFitness(){
    _generation->sortByFitness();
}

void CrystalSearch::_writeXYZ(double myrCut, const char* filename){
    _setUp(&(_generation->individuals[0]),true);
    
    double a1l, a2l, a3l;
    if (myrCut == 0.0) myrCut = 4.0;
    
    a1l = _a1.getLength(); a2l = _a2.getLength(); a3l = _a3.getLength();
    
    const int iMax = ceil(myrCut/a1l) + 1;
    const int jMax = ceil(myrCut/a2l) + 1;
    const int kMax = ceil(myrCut/a3l) + 1;
    
    
    std::ofstream myinfo("box_config.dat");
    //myinfo << "Box info for python script\n";
    myinfo << (2*iMax+1)*_a1.x << "\t" << (2*jMax+1)*_a2.x << "\t" << (2*jMax+1)*_a2.y << "\t"
    << (2*kMax+1)*_a3.x << "\t" << (2*kMax+1)*_a3.y << "\t" << (2*kMax+1)*_a3.z << "\n";
    myinfo.close();
    
    std::ofstream myxyz(filename);
    myxyz << ((2*iMax + 1)*(2*jMax + 1 )*(2*kMax + 1)*_mymolecules.size())*_mymolecules.begin()->first.size() << "\n"; //assumes all molecules are the same size... it might make sense to have a counter here
    myxyz << (2*iMax+1)*_a1.x << "\t" << (2*jMax+1)*_a2.x << "\t" << (2*jMax+1)*_a2.y << "\t"
    << (2*kMax+1)*_a3.x << "\t" << (2*kMax+1)*_a3.y << "\t" << (2*kMax+1)*_a3.z << "\tbox_dim\n";
    
    int ncount=0;
    for (int i=-iMax; i<= iMax; i++){
        for (int j=-jMax; j<= jMax; j++ ){
            for (int k=-kMax; k<=kMax; k++){
                const double3 site_ijk = double3(i*_a1.x + j*_a2.x + k*_a3.x,
                                                 i*_a1.y + j*_a2.y + k*_a3.y,
                                                 i*_a1.z + j*_a2.z + k*_a3.z);
                
                for (const auto& m : _mymolecules){
                    for (auto& mx : m.first){
                        myxyz << m.second << "\t" << mx.x + site_ijk.x << "\t" << mx.y + site_ijk.y << "\t" << mx.z + site_ijk.z << "\n";
                        ncount++;
                    }
                }
                
            }
        }
    }
    std::cout << "Total lines written is " << ncount << "\n";
    
    
}


void CrystalSearch::printBestIndividual(){
    _generation->sortByFitness();
    individual* best_ind = &(_generation->individuals[0]);
    
    std::ofstream bPInfo("bestSolutionParameter.dat",std::ios_base::app | std::ios_base::out);
    std::ofstream bVInfo("bestSolutionVector.dat",std::ios_base::app | std::ios_base::out);
    
    
    //write solution parameter values
    bPInfo << _gencount << "\t";
    for (auto& m : best_ind->genes) bPInfo << m.getValue() << "\t";
    bPInfo << getBittage() << "\n";
    
    //write solution vector values
    _setUp(best_ind);
    
    bVInfo << _gencount << "\t";
    bVInfo << _a1.x << "\t" << _a1.y << "\t" << _a1.z << "\t";
    bVInfo << _a2.x << "\t" << _a2.y << "\t" << _a2.z << "\t";
    bVInfo << _a3.x << "\t" << _a3.y << "\t" << _a3.z << "\t";
    
    for (auto& s : _dcm) bVInfo << s.x << "\t" << s.y << "\t" << s.z << "\t";
    for (auto& s : _axangle) bVInfo << s.coords.x << "\t" << s.coords.y << "\t" << s.coords.z << "\t" << s.scalar << "\t";
    
    if (_myconfig->s_flag) bVInfo << (best_ind->genes.back()).getValue();
    bVInfo << "\n";
    
    bPInfo.close();
    bVInfo.close();
}

void CrystalSearch::_printIndividual(unsigned int i){
    
    _ind = &(_generation->individuals[i]);
    
    //write solution parameter values
    *solPInfo << _gencount << "\t";
    for (auto& m : _ind->genes) *solPInfo << m.getValue() << "\t";
    *solPInfo << "\n";
    
    //write solution vector values
    _setUp(_ind);
    
    *solVInfo << _gencount << "\t";
    *solVInfo << _a1.x << "\t" << _a1.y << "\t" << _a1.z << "\t";
    *solVInfo << _a2.x << "\t" << _a2.y << "\t" << _a2.z << "\t";
    *solVInfo << _a3.x << "\t" << _a3.y << "\t" << _a3.z << "\t";
    
    for (auto& s : _dcm) *solVInfo << s.x << "\t" << s.y << "\t" << s.z << "\t";
    for (auto& s : _axangle) *solVInfo << s.coords.x << "\t" << s.coords.y << "\t" << s.coords.z << "\t" << s.scalar << "\t";
    
    if (_myconfig->s_flag) *solVInfo << (_ind->genes.back()).getValue();
    *solVInfo << "\n";
    
    
}

void CrystalSearch::_setUp(individual* ind, bool rotate_flag){
     _genes = &(ind->genes);
    
    //lattice operations
    _lattice.updateVectors(_genes->at(0).getValue(),
                           _genes->at(1).getValue(),
                           _genes->at(2).getValue(),
                           _genes->at(3).getValue(),
                           _genes->at(4).getValue());
    _lattice.setDensity(_density/_myconfig->molecules.size());
    
    _a1 = _lattice.getLatticeVector(0);
    _a2 = _lattice.getLatticeVector(1);
    _a3 = _lattice.getLatticeVector(2);
    
    unsigned int k=5; //remove hard coding
    //get other basis point
    for (unsigned int l=1; l<_dcm.size(); l++){
        for (auto& m : _cms) m = (*_genes)[k++].getValue();
        _dcm[l] = double3(_a1.x*_cms[0] + _a2.x*_cms[1] + _a3.x+_cms[2],
                          _a1.y*_cms[0] + _a2.y*_cms[1] + _a3.y*_cms[2],
                          _a1.z*_cms[0] + _a2.z*_cms[1] + _a3.z*_cms[2]);
    }
        
    //get orientation for other base pts.
    for (unsigned int l=0; l<_quat.size(); l++){
        for (unsigned int i=0; i < 4; i++ ) _qts[i] = (*_genes)[k++].getValue();
        normalize(_qts[0], _qts[1], _qts[2], _qts[3]);
        _quat[l] = double4(_qts[0],_qts[1],_qts[2],_qts[3]);
        _axangle[l] = axisangle(_quat[l]);
        
    }
    
    //update molecules (use a flag here)
    //maybe allow 
    if (rotate_flag){
        int combit = (int) (ONEMINUSEPSILON*(_myconfig->ncombinations)*(_genes->back().getValue()));
        //test the distribution of sampling, cheat to ensure range is correctly sampled
        assert(combit < _myconfig->combinations.size());
        assert(_quat.size() >= _myconfig->combinations[combit].size());
        _mymolecules.resize(_myconfig->combinations[combit].size()); //maybe do this differently
        unsigned int l=0;
        for (const auto& str : _myconfig->combinations[combit]){
            _mymolecules[l].second = _myconfig->molecules[_typemap[str]].second;
            _mymolecules[l].first = rotate(_quat[l],_myconfig->molecules[_typemap[str]].first);
            zeroCOM(_mymolecules[l].first); //fix at zero
            for (auto&t : _mymolecules[l].first) t = t + _dcm[l]; //translate
            l++;

        }
        
    }
    
    
    
    
    
}

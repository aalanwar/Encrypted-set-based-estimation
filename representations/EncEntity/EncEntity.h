#pragma once

#include <vector>
#include <typedef.h>
#include <random>
#include <paillier.hh>
#include <EncZonotope/EncZonotope.h>
#include <EncStrip/EncStrip.h>



namespace EncEst {

class EncEntity {
  public:
    vector<EncStrip> mEncStripList;
    std::default_random_engine mGen;
  	std::uniform_real_distribution<double> mDist;
    Paillier* mP;
    EncZonotope mEncZ0;


    EncEntity(const matrix_t<double>& HH, const matrix_t<double>& RR, Paillier* p);
    //void measure(const vector_t<double>& x);


};

}  // namespace

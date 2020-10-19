

#pragma once

#include "../typedef.h"
#include <paillier.hh>
#include <vector>
#include <EncStrip/EncStrip.h>
#include <mat.h>
#include <eigen3/Eigen/Dense>
#include <map>
#include <iostream>
namespace EncEst {

class EncConZonotope {
  public:
		size_t mDimension;
		vector_t<encnum> mEncCenter;
		matrix_t<double> mGenerators;
		matrix_t<double> mA;
		vector_t<encnum> mEncb;
		Paillier* mP;

		
	EncConZonotope(const vector_t<encnum>& enccenter, const matrix_t<double>& generators, const matrix_t<double>& A,const vector_t<encnum>& encb ,Paillier* p);
	
    EncConZonotope intersect(const EncConZonotope& ey);
    
    EncConZonotope intersect(const std::vector<EncStrip>& encstrips);

		//overloaded operators
	EncConZonotope operator+(const EncConZonotope& ez);
	static EncConZonotope intersectMany(const std::vector<EncConZonotope>& conzonol) ;

};
//global overloaded operators
//EncZonotope operator*(const matrix_t<double>& mat, EncZonotope &ez);

}  // namespace

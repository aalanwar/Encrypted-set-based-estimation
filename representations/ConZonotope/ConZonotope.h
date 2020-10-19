

#pragma once

#include <typedef.h>
#include <vector>
#include <Strip/Strip.h>
#include <mat.h>
#include <eigen3/Eigen/Dense>
#include <map>
#include <iostream>
namespace EncEst {

class ConZonotope {
  public:
		size_t mDimension;
		vector_t<double> mCenter;
		matrix_t<double> mGenerators;
		matrix_t<double> mA;
		vector_t<double> mb;
		ConZonotope(const vector_t<double>& center, const matrix_t<double>& generators,const matrix_t<double>& A,const vector_t<double>& b);

    ConZonotope intersect(const ConZonotope& ey);
    
    ConZonotope intersect(const std::vector<Strip>& strips);

		//overloaded operators
	ConZonotope operator+(const ConZonotope& ez);

};
}  // namespace

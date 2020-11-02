#pragma once

#include <vector>
#include <mat.h>
#include <eigen3/Eigen/Dense>
#include <map>
#include <iostream>

#include "Strip/Strip.h"
#include "typedef.h"

namespace EncEst 
{
	class ConZonotope 
	{
	public:
		size_t mDimension;
		vector_t<double> mCenter;
		matrix_t<double> mGenerators;
		matrix_t<double> mA;
		vector_t<double> mb;

		ConZonotope(const vector_t<double>& center, const matrix_t<double>& generators,const matrix_t<double>& A,const vector_t<double>& b);

		ConZonotope intersect(const ConZonotope& ey);
		
		ConZonotope intersect(const std::vector<Strip>& strips);

		ConZonotope operator+(const ConZonotope& ez);
	};
}

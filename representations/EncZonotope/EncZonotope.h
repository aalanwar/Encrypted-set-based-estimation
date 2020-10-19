

#pragma once

#include "../typedef.h"
#include <paillier.hh>
#include <vector>
#include <EncStrip/EncStrip.h>
#include <mat.h>
#include <eigen3/Eigen/Dense>
#include <map>
#include <iostream>

#include <queue>
#include <chrono>
#include <limits>
#include <algorithm>
#include<cmath>
 
using namespace std::chrono;
 
typedef std::pair<double, int> idx_pair;

namespace EncEst {

class EncZonotope {
  public:
		size_t mDimension;
		vector_t<encnum> mEncCenter;
		matrix_t<double> mGenerators;
		Paillier* mP;

		EncZonotope(const vector_t<encnum>& enccenter, const matrix_t<double>& generators, Paillier* p);
		EncZonotope(const EncZonotope& rhs);
    EncZonotope intersect(const EncZonotope& ez);
    static EncZonotope intersectMany(const std::vector<EncZonotope>& ezs);
    EncZonotope intersect(const std::vector<EncStrip>& encstrips);
	EncZonotope reduce(const int order);
	std::vector<idx_pair> getNSmallest(std::vector<double> const& data, int n);
		//overloaded operators
		EncZonotope operator+(const EncZonotope& ez);
		//EncZonotope operator=(const EncZonotope& rhs);
	bool equal(double value1, double value2);
    std::vector<matrix_t<double>> computeLambdas(const std::vector<EncStrip>& encstrips);
    static vector_t<double> computeWeights(const std::vector<EncZonotope>& enczonos);
};
//global overloaded operators
EncZonotope operator*(const matrix_t<double>& mat, EncZonotope &ez);


}  // namespace

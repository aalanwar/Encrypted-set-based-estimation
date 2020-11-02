#pragma once

#include <vector>
#include <mat.h>
#include <eigen3/Eigen/Dense>
#include <map>
#include <iostream>
#include <queue>
#include <chrono>
#include <limits>
#include <algorithm>
#include <cmath>

#include "typedef.h"
#include "paillier.hh"
#include "EncStrip/EncStrip.h"

using namespace std::chrono;
 
typedef std::pair<double, int> idx_pair;

namespace EncEst 
{
	/// @brief An encrypted Zonotope respresentation.
	class EncZonotope 
	{
	public:
		size_t mDimension;
		vector_t<encnum> mEncCenter;
		matrix_t<double> mGenerators;
		Paillier* mP;

	public:
		/// @brief Constructs a new encrypted Zonotope.
		/// @param enccenter The center vector of the encrypted Zonotope.
		/// @param generators The generator matrix of the encrypted Zonotope.
		/// @param p The Paillier encryption scheme.
		EncZonotope(const vector_t<encnum>& enccenter, const matrix_t<double>& generators, Paillier* p);
		EncZonotope(const EncZonotope& rhs);
		
		/// @brief Calculates the intersection of the encrypted Zonotope with another encrypted Zonotope.
		/// @param ez A collection of encrypted Strips to intersect the encrypted Zonotope with.
		/// @return A new encrypted Zonotope which is a tight overapproximation of the intersection between the original two Zonotopes.
		EncZonotope intersect(const EncZonotope& ez);
		
		/// @brief Calculates the intersection of the encrypted Zonotope with a collection of encrypted Strips.
		/// @param strips A collection of encrypted Strips to intersect the encrypted Zonotope with.
		/// @return A new encrypted Zonotope which is a tight overapproximation of the intersection between the original encrypted Zonotope and the encrypted Strips.
		EncZonotope intersect(const std::vector<EncStrip>& encstrips);
		
		/// @brief Performs a reduce operation on the encrypted Zonotope.
		/// @param order TODO
		/// @return A new encrypted Zonotope which has smaller generator matrix.
		EncZonotope reduce(const int order);
		
		
		/// @brief TODO
		/// @param data TODO
		/// @param n TODO
		/// @return TODO
		std::vector<idx_pair> getNSmallest(std::vector<double> const& data, int n);

		/// @brief Adds two encrypted Zonotopes.
		/// @param ez Another encrypted Zonotope.
		/// @return A new encrypted Zonotope which is equal to the Minkowsky sum of the two encrypted Zonotopes. 
		EncZonotope operator+(const EncZonotope& ez);

		//EncZonotope operator=(const EncZonotope& rhs);
		bool equal(double value1, double value2);
		virtual std::vector<matrix_t<double>> computeLambdas(const std::vector<EncStrip>& encstrips);

	public:
		/// @brief TODO
		/// @param data TODO
		/// @param n TODO
		/// @return TODO
		static EncZonotope intersectMany(const std::vector<EncZonotope>& ezs);

		/// @brief TODO
		/// @param data TODO
		/// @param n TODO
		/// @return TODO
		static vector_t<double> computeWeights(const std::vector<EncZonotope>& enczonos);
	};

	/// @brief Multiplies a Matrix with an encrypted Zonotope
	/// @param mat A matrix with valid dimensions.
	/// @param ez The encrypted Zonotopoe.
	/// @return A new encrypted Zonotope, which is the TODO...
	EncZonotope operator*(const matrix_t<double>& mat, EncZonotope &ez);

}  // namespace

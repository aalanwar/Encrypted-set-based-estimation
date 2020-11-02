#pragma once

#include "typedef.h"
#include "Strip/Strip.h"

namespace EncEst 
{
	/// @brief A compact representation of a high-dimensional set, that are closed under Linear-Transformations and Minkowski addition.
	class Zonotope 
	{
	public:
		size_t mDimension;
		vector_t<double> mCenter;
		matrix_t<double> mGenerators;

		/// @brief Constructs a new Zonotope.
		/// @param center The center vector of the Zonotope.
		/// @param generators The generator matrix of the Zonotope. 
		Zonotope(const vector_t<double>& center, const matrix_t<double>& generators);

		/// @brief TODO
		virtual std::vector<matrix_t<double>> computeLambdas(const std::vector<Strip>& strips);
		
		/// @brief Calculates the intersection of the Zonotope with a collection of Strips.
		/// @param strips A collection of Strips to intersect the Zonotope with.
		/// @return A new Zonotope which is a tight overapproximation of the intersection between the original Zonotope and the Strips.
		Zonotope intersect(const std::vector<Strip>& strips);
	};
}

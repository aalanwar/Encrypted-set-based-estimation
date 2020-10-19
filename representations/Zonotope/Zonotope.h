

#pragma once

#include <typedef.h>
#include <Strip/Strip.h>
namespace EncEst {

class Zonotope {
  public:
		size_t mDimension;
		vector_t<double> mCenter;
		matrix_t<double> mGenerators;

		Zonotope(const vector_t<double>& center, const matrix_t<double>& generators);
		std::vector<matrix_t<double>> computeLambdas(const std::vector<Strip>& strips);
		Zonotope intersect(const std::vector<Strip>& strips);
};
}  // namespace

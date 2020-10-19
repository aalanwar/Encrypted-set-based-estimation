#pragma once
#include <eigen3/Eigen/Dense>

namespace EncEst {
	template <typename Number>
	using vector_t = Eigen::Matrix<Number, Eigen::Dynamic, 1>;
	
	template <typename Number>

	using matrix_t = Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>;
}

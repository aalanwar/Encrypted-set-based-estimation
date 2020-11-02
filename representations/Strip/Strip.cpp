#include <Strip/Strip.h>
#include <stdexcept>

namespace EncEst 
{
	Strip::Strip(const matrix_t<double>& H, const vector_t<double>& y, const vector_t<double> r)
		: mH(H)
		, mY(y)
		, mR(r)
	{
		// TODO: there are propably better ways to check for good dimensions (look into Eigen source...)
		// assert(H.rows() == 1);
		// assert(y.size() == 1);
		// assert(r.size() == 1);

		// TODO: assert r,y € 1x1 = scalar.
		// For a strip it holds that r,y € 1x1 and H € 1xn with n = dim(zonotope)
	}

	Strip::Strip(const matrix_t<double>& H, const double y, const double r)
		: mH{H}
		, mY{1,1}
		, mR{1,1}
	{
		// TODO: makes no sense to check dimensionality here because mH, mY, mR are public so, can be changed anyways ...
		// assert(H.rows() == 1);

		mY << y;
		mR << r;
	}
}

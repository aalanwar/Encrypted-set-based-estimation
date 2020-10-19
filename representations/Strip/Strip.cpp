#include <Strip/Strip.h>
#include <stdexcept>
namespace EncEst {

Strip::Strip(const matrix_t<double>& H, const vector_t<double>& y, const vector_t<double> r):
 														mH(H),
														mY(y),
														mR(r)
														{}

}//end namespace

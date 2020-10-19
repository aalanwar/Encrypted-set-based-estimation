#pragma once

#include <vector>
#include <typedef.h>

namespace EncEst {

class Strip {
  public:
    matrix_t<double> mH;
    vector_t<double> mR;
    vector_t<double> mY;

    Strip(const matrix_t<double>& H, const vector_t<double>& y, const vector_t<double> r);
};

}  // namespace

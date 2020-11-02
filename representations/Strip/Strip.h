#pragma once

#include <vector>
#include "typedef.h"

namespace EncEst 
{
  /// @brief A compact representation of a measurement state.
  class Strip 
  {
  public:
      matrix_t<double> mH;
      vector_t<double> mY;
      vector_t<double> mR;

      /// @brief Construct a new Stip.
      /// @param H Measurement matrix.
      /// @param y Measurement of a sensor.
      /// @param r Measurement uncertainty vector for all sensors.
      Strip(const matrix_t<double>& H, const vector_t<double>& y, const vector_t<double> r);
      
      // SJ: TODO
      Strip(const matrix_t<double>& H, const double y, const double r);
  };
}

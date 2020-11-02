#pragma once

#include "typedef.h"
#include "paillier.hh"
#include "Strip/Strip.h"

namespace EncEst 
{
  class EncStrip 
  {
    public:
      Paillier* mP;
      matrix_t<double> mH;
      vector_t<double> mR;
      vector_t<encnum> mEncY;

      EncStrip(const matrix_t<double>& H, const vector_t<encnum>& ency, const vector_t<double> r, Paillier* p);
      EncStrip(const Strip& strip, Paillier* p);

      //std::vector<Point<double>> vertices( const matrix_t<double>& = matrix_t<double>::Zero(0,0) ) const;
  };
}

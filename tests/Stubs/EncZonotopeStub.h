#pragma once

#include "EncZonotope/EncZonotope.h"

namespace EncEst
{
    namespace Test
    {
        class EncZonotopeStub : public EncZonotope
        {
        public:
            // TODO: Set lambdas as public member instead of calculating them ...
            std::vector<matrix_t<double>> lambdas {};

            EncZonotopeStub(const vector_t<encnum>& enccenter, const matrix_t<double>& generators, Paillier* p)
                : EncZonotope(enccenter, generators, p)
            {}

            std::vector<matrix_t<double>> computeLambdas(const std::vector<EncStrip>& encstrips) override
            {
                if(lambdas.size() > 0)
                    return lambdas;

                assert(encstrips.size() > 0);

                std::vector<matrix_t<double>> v(encstrips.size());

                int n = mDimension;
                int p = encstrips[0].mEncY.rows();

                for(size_t i = 0; i < encstrips.size(); ++i)
                    v[i] = matrix_t<double>::Ones(n, p);

                return v;
            }
        };
    }
}
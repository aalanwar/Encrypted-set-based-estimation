#pragma once

#include "Zonotope/Zonotope.h"

namespace EncEst
{
    namespace Test
    {
        class ZonotopeStub : public Zonotope
        {
        public:
            // TODO: Set lambdas as public member instead of calculating them ...
            std::vector<matrix_t<double>> lambdas {};

            ZonotopeStub(const vector_t<double>& center, const matrix_t<double>& generators)
                : Zonotope(center, generators)
            {}

            std::vector<matrix_t<double>> computeLambdas(const std::vector<Strip>& strips) override
            {
                if(lambdas.size() > 0)
                    return lambdas;

                assert(strips.size() > 0);

                std::vector<matrix_t<double>> v(strips.size());

                int n = mDimension;
                int p = strips[0].mY.rows();

                for(size_t i = 0; i < strips.size(); ++i)
                    v[i] = matrix_t<double>::Ones(n, p);

                return v;
            }
        };
    }
}
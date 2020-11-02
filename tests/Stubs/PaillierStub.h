#pragma once

#include "paillier.hh"

namespace EncEst
{
    namespace Test
    {
        class PaillierStub : public Paillier
        {
        public:
            const static inline size_t DEFAULT_SIZE_RAND_VALUES {10};
            std::vector<mpz_class> FixedRandomValues{DEFAULT_SIZE_RAND_VALUES};

            PaillierStub(const std::vector<mpz_class> &pk, gmp_randstate_t state)
                : Paillier(pk, state)
            {}

            void InitFixedRandomValues(const size_t size = DEFAULT_SIZE_RAND_VALUES)
            {
                size_t count = size - FixedRandomValues.size();
                mpz_class r;

                for (size_t i = 0; i < count; ++i) 
                {
                    mpz_urandomm(r.get_mpz_t(),_randstate,n.get_mpz_t());
                    FixedRandomValues.push_back(mpz_class_powm(g,n*r,n2));
                }
            }

            void ResetRandomList(const size_t size = DEFAULT_SIZE_RAND_VALUES)
            {
                if(FixedRandomValues.size() < size)
                    InitFixedRandomValues(size);

                this->rqueue.clear();

                for(size_t i = 0; i < size; ++i)
                    this->rqueue.push_back(FixedRandomValues[i]);
            }
        
        };
    }
}
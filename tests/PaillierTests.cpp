#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <eigen3/Eigen/Dense>

#include "TestUtils.h"
#include <gmp.h>
#include <gmpxx.h>

namespace EncEst
{
    namespace Test
    {
        using namespace testing;
        using namespace EncEst::Test::Utils;

        // Fixture
        class TPaillier : public Test
        {
        public:
            PaillierUtil paillierUtil {};

            void SetUp() override 
            {}
        };

        // Tests
        TEST_F(TPaillier, EqualDoubleEncryptsToSameValue)
        {
            paillierUtil.p->ResetRandomList();
            encnum en1 = paillierUtil.p->encrypt_f(100);

            paillierUtil.p->ResetRandomList();
            encnum en2 = paillierUtil.p->encrypt_f(100);

            ASSERT_THAT(en1, Eq(en2));
        }


        TEST_F(TPaillier, UnequalDoubleEncryptsToDifferentValue)
        {
            paillierUtil.p->ResetRandomList();
            encnum en1 = paillierUtil.p->encrypt_f(100);

            paillierUtil.p->ResetRandomList();
            encnum en2 = paillierUtil.p->encrypt_f(200);

            ASSERT_THAT(en1, Ne(en2));
        }

        // TEST_F(TPaillier, TODO)
        // {
        //     mpz_class c1, c2, c3, c4;
        //     c1.set_str("12345678901234566890123455678901234567890123456689012345567890", 10);
        //     c2 = 120;
        //     c3.set_str("12345678901234566890123455678901234567890123456689012345568010", 10);
        //     c4 = c1 + c2;

        //     std::cout << c1 << std::endl;
        //     std::cout << c2 << std::endl;
        //     std::cout << c3 << std::endl;
        //     std::cout << c4 << std::endl;

        //     ASSERT_THAT(mpz_cmp(c3.get_mpz_t(), c4.get_mpz_t()), Eq(0));
        //     ASSERT_THAT(c3 == c4, Eq(true));
        // }
    }
}
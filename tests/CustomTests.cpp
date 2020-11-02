#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <eigen3/Eigen/Dense>

#include "TestUtils.h"

namespace EncEst
{
    namespace Test
    {
        using namespace testing;
        using namespace EncEst::Test::Utils;

        // Fixture
        class TCustom : public Test
        {
        public:
            void SetUp() override 
            {}
        };

        TEST_F(TCustom, MatrixDefaultsToDouble)
        {
            auto A = Matrix::init<3,3>({1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0});
            auto B = Matrix::init<3,3>({1,1,1,1,1,1,1,1,1});

            matrix_t<double> C(3,3);
            C << 1,1,1,
                 1,1,1,
                 1,1,1;

            auto D = Matrix::init<3,3>({9,9,9,9,9,9,9,9,9}); 

            ASSERT_THAT(A*B*C, Eq(D));
        }

        // Tests
        TEST_F(TCustom, CreatesValid3x3Matrix)
        {
            auto A = Matrix::init<3,3>({1,2,3,4,5,6,7,8,9});

            matrix_t<double> B(3,3);
            B << 1,2,3,
                 4,5,6,
                 7,8,9;

            ASSERT_THAT(A, Eq(B));
        }

        TEST_F(TCustom, CreatesValid5x3Matrix)
        {
            auto A = Matrix::init<5,3>({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15});

            matrix_t<double> B(5,3);
            B <<  1,  2,  3,
                  4,  5,  6,
                  7,  8,  9,
                 10, 11, 12,
                 13, 14, 15;

            ASSERT_THAT(A, Eq(B));
        }

        TEST_F(TCustom, CreatesValid3x5Matrix)
        {
            auto A = Matrix::init<3,5>({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15});

            matrix_t<double> B(3,5);
            B <<  1,  2,  3,  4,  5,
                  6,  7,  8,  9, 10,
                 11, 12, 13, 14, 15;

            ASSERT_THAT(A, Eq(B));
        }

        TEST_F(TCustom, OnlyUsesFirstXElementsForMatrixInitialization)
        {
            auto A = Matrix::init<3,3>({1,2,3,4,5,6,7,8,9,10,11,12});

            matrix_t<double> B(3,3);
            B << 1,2,3,
                 4,5,6,
                 7,8,9;

            ASSERT_THAT(A, Eq(B));
        }

        TEST_F(TCustom, AssertsOnWrongMatrixInitializerSize)
        {
            auto fail = []() {Matrix::init<3,3>({1,2,3,4,5,6,7});};
            ASSERT_DEATH(fail(), "Assertion");
        }

        TEST_F(TCustom, CreatesValidRowVector)
        {
            auto a = Vector::row<3>({1,2,3});

            matrix_t<double> b(1,3);
            b << 1,2,3;

            ASSERT_THAT(a, Eq(b));
        }

        TEST_F(TCustom, CreatesValidColumnVector)
        {
            auto a = Vector::col<3>({1,2,3});

            vector_t<double> b(3,1);
            b << 1,
                 2,
                 3;

            ASSERT_THAT(a, Eq(b));
        }

        TEST_F(TCustom, CreatesValid1x3Vector)
        {
            auto a = Vector::init<1,3>({1,2,3});

            matrix_t<double> b(1,3);
            b << 1,2,3;

            ASSERT_THAT(a, Eq(b));
        }

        TEST_F(TCustom, CreatesValid3x1Vector)
        {
            auto a = Vector::init<3,1>({1,2,3});

            vector_t<double> b(3,1);
            b << 1,2,3;

            ASSERT_THAT(a, Eq(b));
        }

        TEST_F(TCustom, OnlyUsesFirstXElementsForVectorInitialization)
        {
            auto a = Vector::init<1,3>({1,2,3,4,5,6});

            matrix_t<double> b(1,3);
            b << 1,2,3;

            ASSERT_THAT(a, Eq(b));
        }

        TEST_F(TCustom, AssertsOnInvalidVectorTemplateValues)
        {
            auto fail = [](){Vector::init<3,3>({1,2,3,4,5,6,7,8,9});};
            ASSERT_DEATH(fail(), "Assertion");  // gTest expects a regular expression to evaluate...
        }

        TEST_F(TCustom, AssertsOnWrongInitializerSize)
        {
            auto fail = [](){Vector::row<5>({1,2,3});};
            ASSERT_DEATH(fail(), "Assertion");
        }
    }
}
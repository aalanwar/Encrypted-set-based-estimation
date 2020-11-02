#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <eigen3/Eigen/Dense>

#include "Strip/Strip.h"

namespace EncEst
{
    namespace Test
    {
        using namespace testing;

        // Fixture
        class TStrip : public Test
        {
        public:
            vector_t<double> y{1,1};
            vector_t<double> r{1,1};
            matrix_t<double> H{1,3};
            Strip* strip;

            void SetUp() override 
            {
                H << 1,2,3;
                y << 1;
                r << 1;

                strip = new Strip{H,y,r};
            }

            void TearDown()
            {
                delete strip;
            }
        };


        // Test
        TEST_F(TStrip, HasCorrectMeasurementMatrix)
        {
            ASSERT_THAT(strip->mH, Eq(H));
        }

        TEST_F(TStrip, HasCorrectMeassurement)
        {
            ASSERT_THAT(strip->mY, Eq(y));   
        }

        TEST_F(TStrip, HasCorrectUncertainty)
        {
            ASSERT_THAT(strip->mR, Eq(r));   
        }
    
        TEST_F(TStrip, AssertWhenMeasurementMatrixHasInvalidDimension)
        {
            matrix_t<double> Hx {3,2};
            Hx << 1,2,3,4,5,6;
            
            strip->mH = Hx;

            ASSERT_THAT(true, Eq(false));   
        }

        TEST_F(TStrip, AssertWhenMeasurementHasInvalidDimension)
        {
            ASSERT_THAT(true, Eq(false));   
        }

        TEST_F(TStrip, AssertWhenUncertaintyHasInvalidDimension)
        {
            ASSERT_THAT(true, Eq(false));   
        }
    }
}











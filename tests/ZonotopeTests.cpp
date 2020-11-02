#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <eigen3/Eigen/Dense>

#include "Zonotope/Zonotope.h"
#include "Strip/Strip.h"
#include "TestUtils.h"
#include "Stubs/ZonotopeStub.h"

namespace EncEst
{
    namespace Test
    {
        using namespace testing;
        using namespace EncEst::Test::Utils;

        // Fixture
        class TZonotope : public Test
        {
        public:
            Vector3d v1 {Vector3d::Ones()};
            Matrix3d M1 {Matrix3d::Identity()};
            Zonotope zonotope1 {v1, M1};
            std::vector<Strip> strips1 {};

            void SetUp() override 
            {
                vector_t<double> y(1,1);
                y << 1;
                vector_t<double> R(1,1);
                R << 1;
                matrix_t<double> H(1,3);
                H << 1,2,3;

                strips1.assign 
                ({
                    Strip {H * 1, y * 1, R * 1},
                    Strip {H * 2, y * 2, R * 2},
                    Strip {H * 3, y * 3, R * 3},
                });

                std::cout << strips1.size() << std::endl;
            }
        };


        // Tests
        TEST_F(TZonotope, HasCorrectDimension)
        {
            ASSERT_THAT(zonotope1.mDimension, Eq(v1.rows()));
        }

        TEST_F(TZonotope, HasCorrectCenter)
        {
            ASSERT_THAT(zonotope1.mCenter, Eq(v1));
        }

        TEST_F(TZonotope, HasCorrectGenerators)
        {
            ASSERT_THAT(zonotope1.mGenerators, Eq(M1));  
        }
        
        TEST_F(TZonotope, IntersectionWithStripsReturnsNewZonotopeWithCorrectDimension)
        {
            ZonotopeStub z1 {v1, M1};
            Zonotope z2 = z1.intersect(strips1);

            ASSERT_THAT(z2.mDimension, Eq(z1.mDimension));
        }

        TEST_F(TZonotope, IntersectionWithStripsReturnsNewZonotopeWithCorrectCenter)
        {
            ZonotopeStub z1 {v1, M1};
            Zonotope z2 = z1.intersect(strips1);

            vector_t<double> v(3,1);
            v << -29, -29, -29;

            ASSERT_THAT(z2.mCenter, Eq(v));
        }

        TEST_F(TZonotope, IntersectionWithStripsReturnsNewZonotopeWithCorrectGenerators)
        {
            ZonotopeStub z1 {v1, M1};
            Zonotope z2 = z1.intersect(strips1);

            matrix_t<double> M(3,6);
            M << -5, -12, -18, 1, 2, 3,
                 -6, -11, -18, 1, 2, 3,
                 -6, -12, -17, 1, 2, 3;

            ASSERT_THAT(z2.mGenerators, Eq(M));
        }

        TEST_F(TZonotope, AssertsWhenCenterAndGeneratorHaveInvalidDimensions)
        {
            vector_t<double> c(3,1);
            c << 1,2,3;
            
            matrix_t<double> G(2,4);
            G << 1,2,
                 3,4,
                 5,6,
                 7,8;

            ASSERT_DEATH(Zonotope(c,G), Eq("TODO"));
        }

        TEST_F(TZonotope, AssertsWhenIntersectingInvalidDimensionsForZonotopesAndStrips)
        {
            matrix_t<double> H(1,4);
            H << 1,2,3,4;

            std::vector<Strip> strips 
            {
                Strip {H, 1, 1}
            };

            // Zonotope with dim == 3 and Strips with dim == 4
            ASSERT_DEATH(zonotope1.intersect(strips), Eq("TODO"));
        }

        TEST_F(TZonotope, ConstructsCorrectNumberOfLambdas)
        {
            auto lambdas = zonotope1.computeLambdas(strips1);

            ASSERT_THAT(lambdas.size(), Eq(strips1.size()));
        }

        TEST_F(TZonotope, ConstructsLambdasWithCorrectDimensions)
        {
            auto lambdas = zonotope1.computeLambdas(strips1);

            int cols {0};
            for(const auto& l : lambdas)
            {
                cols += l.cols();
                ASSERT_THAT(l.rows(), Eq(strips1[0].mY.rows()));
            }

            ASSERT_THAT(cols, Eq(zonotope1.mDimension));
        }

        TEST_F(TZonotope, ConstructsCorrectLambdas)
        {
            // Check if lambda values are correct.
            ASSERT_THAT(true, Eq(false));
        }
    }
}











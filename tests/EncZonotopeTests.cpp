#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <eigen3/Eigen/Dense>

#include "EncZonotope/EncZonotope.h"
#include "Strip/Strip.h"
#include "TestUtils.h"
#include "Stubs/ZonotopeStub.h"
#include "Stubs/EncZonotopeStub.h"

namespace EncEst
{
    namespace Test
    {
        using namespace testing;
        using namespace EncEst::Test::Utils;

        // Fixture
        class TEncZonotope : public Test
        {
        public:
            Vector3d v1 = Vector::init<3>({1,1,1});
            Vector3d v2 = v1*2;
            Vector3d v3 = v1*3;
            Matrix3d M1 = Matrix3d::Identity();
            Matrix3d M2 = M1*2;
            Matrix3d M3 = M1*3;
            Zonotope zonotope1 {v1, M1};
            Zonotope zonotope2 {v2, M2};
            Zonotope zonotope3 {v3, M3};
            vector_t<double> r1 = Vector::init<1>({1});
            vector_t<double> y1 = Vector::init<1>({1});
            matrix_t<double> H1 = Matrix::init<1,3>({1,2,3});
            std::vector<Strip> strips1 
            {
                Strip {H1 * 1, y1 * 1, r1 * 1},
                Strip {H1 * 2, y1 * 2, r1 * 2},
                Strip {H1 * 3, y1 * 3, r1 * 3}
            };

            PaillierUtil paillierUtil {};

            void SetUp() override 
            {}
        };


        // Tests
        TEST_F(TEncZonotope, HasCorrectDimension)
        {
            auto ez = paillierUtil.p->encrypt(zonotope1);

            ASSERT_THAT(ez.mDimension, Eq(zonotope1.mDimension));
        }

        TEST_F(TEncZonotope, HasCorrectCenter)
        {
            paillierUtil.p->ResetRandomList();
            auto ez = paillierUtil.p->encrypt(zonotope1);

            paillierUtil.p->ResetRandomList();
            auto ev = paillierUtil.p->encryptVector_f(zonotope1.mCenter);

            ASSERT_THAT(ez.mEncCenter, Eq(ev));
        }

        TEST_F(TEncZonotope, HasCorrectGenerators)
        {
            auto ez = paillierUtil.p->encrypt(zonotope1);

            ASSERT_THAT(ez.mGenerators, Eq(zonotope1.mGenerators));  
        }
        
        
        TEST_F(TEncZonotope, CreatesValidCopyOfAnotherEncZonotope)
        {
            EncZonotope ez1 = paillierUtil.p->encrypt(zonotope1);
            EncZonotope ez2 {ez1};

            ASSERT_THAT(ez1.mDimension, Eq(ez2.mDimension));
            ASSERT_THAT(&ez1.mDimension, Ne(&ez2.mDimension));

            ASSERT_THAT(ez1.mEncCenter, Eq(ez2.mEncCenter));
            ASSERT_THAT(&ez1.mEncCenter, Ne(&ez2.mEncCenter));

            ASSERT_THAT(ez1.mGenerators, Eq(ez2.mGenerators));
            ASSERT_THAT(&ez1.mGenerators, Ne(&ez2.mGenerators));
        }



        TEST_F(TEncZonotope, IntersectionWithAnotherEncZonotopeReturnsNewEncZonotopeWithCorrectDimension)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = ez1.intersect(ez2);
            auto z = paillierUtil.pp->decrypt(ez3);

            int d = z.mDimension;
            ASSERT_TRUE(ez1.mDimension == d && ez2.mDimension == d && ez3.mDimension == d);
        }

        TEST_F(TEncZonotope, IntersectionWithAnotherEncZonotopeReturnsNewEncZonotopeWithCorrectCenter)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = ez1.intersect(ez2);
            auto z = paillierUtil.pp->decrypt(ez3);

            auto v = Vector::init<3>({1.2,1.2,1.2});

            std::cout << z << std::endl;
            std::cout << v << std::endl;

            ASSERT_THAT(z.mCenter, Eq(v));
        }

        TEST_F(TEncZonotope, IntersectionWithAnotherEncZonotopeReturnsNewEncZonotopeWithCorrectGenerators)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = ez1.intersect(ez2);
            auto z = paillierUtil.pp->decrypt(ez3);

            auto M = Matrix::init<3,6>
            ({
                0.8,   0,   0, 0.4,   0,   0, 
                  0, 0.8,   0,   0, 0.4,   0, 
                  0,   0, 0.8,   0,   0, 0.4
            });

            ASSERT_THAT(z.mGenerators, Eq(M));
        }

        TEST_F(TEncZonotope, IntersectionCommutativity)
        {
            // SJ: TODO: The problem is the ordering in the generator 
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);

            auto o1 = ez1.intersect(ez2);
            auto o2 = ez2.intersect(ez1);
            
            auto z1 = paillierUtil.pp->decrypt(o1);
            auto z2 = paillierUtil.pp->decrypt(o2);

            std::cout << z1 << std::endl;
            std::cout << z2 << std::endl;

            ASSERT_THAT(z1.mDimension, Eq(z2.mDimension));
            ASSERT_THAT(z1.mCenter, Eq(z2.mCenter));
            ASSERT_THAT(z1.mGenerators, Eq(z2.mGenerators));
        }

        // TODO: Also Test Assoziativity and Commutativity for regular Zonotopes.
        TEST_F(TEncZonotope, IntersectionAssoziativity)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = paillierUtil.p->encrypt(zonotope3);
            auto ez4_o1 = ez1.intersectMany({ez2, ez3});
            auto ez4_o2 = ez2.intersectMany({ez3, ez1});
            auto z1 = paillierUtil.pp->decrypt(ez4_o1);
            auto z2 = paillierUtil.pp->decrypt(ez4_o2);

            std::cout << z1 << std::endl;
            std::cout << z2 << std::endl;

            ASSERT_THAT(z1.mDimension, Eq(z2.mDimension));
            ASSERT_THAT(z1.mCenter, Eq(z2.mCenter));
            ASSERT_THAT(z1.mGenerators, Eq(z2.mGenerators));
        }

        TEST_F(TEncZonotope, IntersectionWithMultipleEncZonotopeReturnsNewEncZonotopeWithCorrectDimension)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = paillierUtil.p->encrypt(zonotope3);
            auto ez4 = ez1.intersectMany({ez2, ez3});
            auto z = paillierUtil.pp->decrypt(ez4);

            int d = z.mDimension;
            ASSERT_TRUE(ez1.mDimension == d && ez2.mDimension == d && ez3.mDimension == d);
        }

        TEST_F(TEncZonotope, IntersectionWithMultipleEncZonotopeReturnsNewEncZonotopeWithCorrectCenter)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = paillierUtil.p->encrypt(zonotope3);
            auto ez4 = ez1.intersectMany({ez2, ez3});
            auto z = paillierUtil.pp->decrypt(ez4);

            std::cout << z << std::endl;

            auto v = Vector::init<3>({1.2,1.2,1.2});

            ASSERT_THAT(z.mCenter, Eq(v));
        }

        TEST_F(TEncZonotope, IntersectionWithMultipleEncZonotopeReturnsNewEncZonotopeWithCorrectGenerators)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = paillierUtil.p->encrypt(zonotope3);
            auto ez4 = ez1.intersectMany({ez2, ez3});  // TODO: use stub and change computeWeights
            auto z = paillierUtil.pp->decrypt(ez4);

            auto M = Matrix::init<3,6>
            ({
                0.8,   0,   0, 0.4,   0,   0, 
                  0, 0.8,   0,   0, 0.4,   0, 
                  0,   0, 0.8,   0,   0, 0.4,
            });

            ASSERT_THAT(z.mGenerators, Eq(M));
        }



        TEST_F(TEncZonotope, IntersectionWithEncStripsReturnsNewEncZonotopeWithCorrectDimension)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto es1 = paillierUtil.p->encrypt(strips1);
            auto ez2 = stub(ez1).intersect(es1);
            auto z = paillierUtil.pp->decrypt(ez2);

            int d = z.mDimension;
            ASSERT_TRUE(ez1.mDimension == d && ez2.mDimension == d);
        }

        TEST_F(TEncZonotope, IntersectionWithEncStripsReturnsNewEncZonotopeWithCorrectCenter)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto es1 = paillierUtil.p->encrypt(strips1);
            auto ez2 = stub(ez1).intersect(es1);
            auto z = paillierUtil.pp->decrypt(ez2);

            auto v = Vector::init<3>({-29, -29, -29});

            std::cout << (z.mCenter - v) << std::endl; 

            ASSERT_TRUE(z.mCenter.isApprox(v));
        }

        TEST_F(TEncZonotope, IntersectionWithEncStripsReturnsNewEncZonotopeWithCorrectGenerators)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto es1 = paillierUtil.p->encrypt(strips1);
            auto ez2 = stub(ez1).intersect(es1);
            auto z = paillierUtil.pp->decrypt(ez2);

            matrix_t<double> M(3,6);
            M << -5, -12, -18, 1, 2, 3,
                 -6, -11, -18, 1, 2, 3,
                 -6, -12, -17, 1, 2, 3;

            ASSERT_THAT(z.mGenerators, Eq(M));
        }


        TEST_F(TEncZonotope, AdditionOfTwoEncZonotopesHasCorrectDimension)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = ez1 + ez2;
            auto z = paillierUtil.pp->decrypt(ez3);

            int d = z.mDimension;
            ASSERT_TRUE(ez1.mDimension == d && ez2.mDimension == d && ez3.mDimension == d);
        }

        TEST_F(TEncZonotope, AdditionOfTwoEncZonotopesHasCorrectCenter)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = ez1 + ez2;
            auto z = paillierUtil.pp->decrypt(ez3);

            ASSERT_THAT(z.mCenter, v1 + v2);
        }

        TEST_F(TEncZonotope, AdditionOfTwoEncZonotopesHasCorrectGenerator)
        {
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = paillierUtil.p->encrypt(zonotope2);
            auto ez3 = ez1 + ez2;
            auto z = paillierUtil.pp->decrypt(ez3);

            std::cout << zonotope1 << std::endl;
            std::cout << zonotope2 << std::endl;
            std::cout << z << std::endl;

            matrix_t<double> M(M1.rows(), M1.cols() + M2.cols());
            M << M1, M2;

            ASSERT_THAT(z.mGenerators, M);
        }


        TEST_F(TEncZonotope, MultiplicationOfEncZonotopeWithMatrixHasCorrectDimension)
        {
            auto M = matrix_t<double>::Identity(3,3) * 0.5;
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = M * ez1;
            auto z = paillierUtil.pp->decrypt(ez2);

            int d = z.mDimension;
            ASSERT_TRUE(ez1.mDimension == d && ez2.mDimension == d);
        }

        TEST_F(TEncZonotope, MultiplicationOfEncZonotopeWithMatrixHasCorrectCenter)
        {
            auto M = matrix_t<double>::Identity(3,3) * 0.5;
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = M * ez1;
            auto z = paillierUtil.pp->decrypt(ez2);

            ASSERT_THAT(z.mCenter, M * zonotope1.mCenter);
        }

        TEST_F(TEncZonotope, MultiplicationOfEncZonotopeWithMatrixHasCorrectGenerator)
        {
            auto M = matrix_t<double>::Identity(3,3) * 0.5;
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = M * ez1;
            auto z = paillierUtil.pp->decrypt(ez2);

            ASSERT_THAT(z.mGenerators, M * zonotope1.mGenerators);
        }


        TEST_F(TEncZonotope, ReductionTests)
        {
            auto M = matrix_t<double>::Identity(3,3) * 0.5;
            auto ez1 = paillierUtil.p->encrypt(zonotope1);
            auto ez2 = M * ez1;
            auto z = paillierUtil.pp->decrypt(ez2);

            ASSERT_THAT(z.mGenerators, M * zonotope1.mGenerators);
        }


        // TEST_F(TEncZonotope, AssertsWhenCenterAndGeneratorHaveInvalidDimensions)
        // {
        //     vector_t<double> c(3,1);
        //     c << 1,2,3;
            
        //     matrix_t<double> G(2,4);
        //     G << 1,2,
        //          3,4,
        //          5,6,
        //          7,8;

        //     ASSERT_DEATH(Zonotope(c,G), Eq("TODO"));
        // }

        // TEST_F(TEncZonotope, AssertsWhenIntersectingInvalidDimensionsForZonotopesAndStrips)
        // {
        //     matrix_t<double> H(1,4);
        //     H << 1,2,3,4;

        //     std::vector<Strip> strips 
        //     {
        //         Strip {H, 1, 1}
        //     };

        //     // Zonotope with dim == 3 and Strips with dim == 4
        //     ASSERT_DEATH(zonotope1.intersect(strips), Eq("TODO"));
        // }

        // TEST_F(TEncZonotope, ConstructsCorrectNumberOfLambdas)
        // {
        //     auto lambdas = zonotope1.computeLambdas(strips1);

        //     ASSERT_THAT(lambdas.size(), Eq(strips1.size()));
        // }

        // TEST_F(TEncZonotope, ConstructsLambdasWithCorrectDimensions)
        // {
        //     auto lambdas = zonotope1.computeLambdas(strips1);

        //     int cols {0};
        //     for(const auto& l : lambdas)
        //     {
        //         cols += l.cols();
        //         ASSERT_THAT(l.rows(), Eq(strips1[0].mY.rows()));
        //     }

        //     ASSERT_THAT(cols, Eq(zonotope1.mDimension));
        // }

        // TEST_F(TEncZonotope, ConstructsCorrectLambdas)
        // {
        //     // Check if lambda values are correct.
        //     ASSERT_THAT(true, Eq(false));
        // }
    }
}
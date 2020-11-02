#pragma once

#include <gm.hh>
#include <NTL/ZZ.h>
#include <memory>
#include <eigen3/Eigen/Dense>
#include <initializer_list>

#include "crypto/paillier.hh"
#include "Zonotope/Zonotope.h"
#include "Strip/Strip.h"
#include "Stubs/PaillierStub.h"
#include "Stubs/ZonotopeStub.h"
#include "Stubs/EncZonotopeStub.h"

namespace EncEst
{
    namespace Test
    {
        namespace Utils
        {
            using namespace NTL;
            using namespace Eigen;

            // TOOD: Remove
            class Values
            {
            public:
                static inline Vector3d v31_1 {Vector3d::Ones() * 1};
                static inline Vector3d v31_2 {Vector3d::Ones() * 2};
                static inline Matrix3d m33_1 {Matrix3d::Identity() * 1};
                static inline Matrix3d m33_2 {Matrix3d::Identity() * 2};

                static inline Zonotope zonotope1 {v31_1, m33_1};
                static inline Zonotope zonotope2 {v31_2, m33_2};
                static inline Strip strip1 {m33_1, v31_1, v31_2};
                static inline Strip strip2 {m33_2, v31_1, v31_2};

            public:
                Values()
                {}

                // static inline Matrix3d m33_1;
                // static inline Zonotope zonotope1 {Eigen::Vector3d::Ones() * 1, Eigen::Matrix3d::Identity()};
                // static inline Zonotope zonotope2 {Eigen::Vector3d::Ones() * 2, Eigen::Matrix3d::Identity()};
            };

            class PaillierUtil
            {
            public:
                shared_ptr<PaillierStub> p;
                shared_ptr<Paillier_priv> pp;

                PaillierUtil(long seed = 999)
                {
                    SetSeed(to_ZZ(seed));
                    
                    gmp_randstate_t randstate;
                    gmp_randinit_default(randstate);
                    gmp_randseed_ui(randstate, seed + 1);

                    auto sk = Paillier_priv::keygen(randstate, 1000, 256);
                    this->pp = std::make_shared<Paillier_priv>(Paillier_priv(sk, randstate));

                    auto pk = this->pp->pubkey();
                    mpz_class nn = pk[0];
                    this->p = make_shared<PaillierStub>(PaillierStub(pk, randstate));
                }
            };
        
            static inline std::ostream& operator<<(std::ostream& stream, const VectorXd& vector)
            {
                for(size_t i = 0; i < vector.rows(); ++i)
                {
                    stream << vector[i] << " ";
                }

                return stream;
            }

            static inline std::ostream& operator<<(std::ostream& stream, const Eigen::MatrixXd matrix)
            {
                for(size_t r = 0; r < matrix.rows(); ++r)
                {
                    for(size_t c = 0; c < matrix.cols(); ++c)
                    {
                        stream << matrix(r,c) << " ";
                    }

                    stream << std::endl;
                }

                return stream;
            }
        
            static inline std::ostream& operator<<(std::ostream& stream, const Strip& strip)
            {
                stream << "[Strip]" << std::endl;
                stream << "measurement sensor:      " << strip.mY << std::endl;
                stream << "measurement uncertainty: " << strip.mR << std::endl;
                stream << "measurement matrix:      " << std::endl << strip.mH;
                
                return stream;
            }

            static inline std::ostream& operator<<(std::ostream& stream, const Zonotope& zonotope)
            {
                stream << "[Zonotope]" << std::endl;
                stream << "dimension:  " << zonotope.mDimension << std::endl;
                stream << "center:     " << zonotope.mCenter << std::endl;
                stream << "generators: " << std::endl << zonotope.mGenerators;

                return stream;
            }
        
            class Matrix2
            {
                public:
                    template<int rows, int cols, typename T = double> 
                    static matrix_t<T> init(const std::initializer_list<T>& args)
                    {
                        assert(args.size() >= rows*cols);

                        matrix_t<T> M(rows, cols);

                        auto it {args.begin()};

                        for(int r = 0; r < rows; ++r)
                        {
                            for(int c = 0; c < cols; ++c)
                            {
                                // SJ: see (http://wiki.ros.org/eigen/Troubleshooting (Syntax for Eigen templates))
                                M.template block<1,1>(r,c) << *it;
                                ++it;
                            }
                        }
                        
                        return M;
                    }
            
                    template<int rows, int cols, typename T = double> 
                    static matrix_t<T> set(const T& d)
                    {
                        matrix_t<T> M(rows, cols);

                        for(int r = 0; r < rows; ++r)
                        {
                            for(int c = 0; c < cols; ++c)
                                M.template block<1,1>(r,c) << d;
                        }
                        
                        return M;
                    }
            };

            class Vector2
            {
                public:
                    template<int rows, int cols, typename T = double>
                    static vector_t<T> init(const std::initializer_list<T>& args)
                    {
                        assert(rows == 1 || cols == 1);
                        assert(args.size() >= rows*cols);

                        vector_t<T> v(rows, cols);

                        auto it {args.begin()};

                        for(int r = 0; r < rows; ++r)
                        {
                            for(int c = 0; c < cols; ++c)
                            {
                                v.template block<1,1>(r,c) << *it;
                                ++it;
                            }
                        }

                        return v;
                    }

                    template<int size, typename T = double>
                    static vector_t<T> row(const std::initializer_list<T>& args)
                    {
                        return Vector2::init<1,size,T>(args);
                    }

                    template<int size, typename T = double>
                    static vector_t<T> col(const std::initializer_list<T>& args)
                    {
                        return Vector2::init<size,1,T>(args);
                    }
            };

            class Vector
            {
                public:
                    // returns matrix_t because vector_t is defined with cols = 1, so we would need to resize vector_t if we want to create a row vector.
                    template<int rows, int cols = 1>
                    static matrix_t<double> init(const std::initializer_list<double>& args)
                    {
                        assert(rows == 1 || cols == 1);
                        assert(args.size() >= rows*cols);

                        matrix_t<double> v(rows, cols);

                        auto it {args.begin()};

                        for(int r = 0; r < rows; ++r)
                        {
                            for(int c = 0; c < cols; ++c)
                            {
                                v.template block<1,1>(r,c) << *it;
                                ++it;
                            }
                        }

                        return v;
                    }

                    template<int rows, int cols = 1>
                    static matrix_t<double> init(const double value)
                    {
                        assert(rows == 1 || cols == 1);
                        return matrix_t<double>::Constant(rows, cols, value);
                    }

                    template<int size>
                    static matrix_t<double> row(const std::initializer_list<double>& args)
                    {
                        return Vector::init<1,size>(args);
                    }

                    template<int size>
                    static matrix_t<double> col(const std::initializer_list<double>& args)
                    {
                        return Vector::init<size,1>(args);
                    }
            };

            class Matrix
            {
                public:
                    template<int rows, int cols> 
                    static matrix_t<double> init(const std::initializer_list<double>& args)
                    {
                        assert(args.size() >= rows*cols);

                        matrix_t<double> M(rows, cols);

                        auto it {args.begin()};

                        for(int r = 0; r < rows; ++r)
                        {
                            for(int c = 0; c < cols; ++c)
                            {
                                // SJ: see (http://wiki.ros.org/eigen/Troubleshooting (Syntax for Eigen templates))
                                M.template block<1,1>(r,c) << *it;
                                ++it;
                            }
                        }
                        
                        return M;
                    }
            };

            static inline EncZonotopeStub stub(const EncZonotope& ez)
            {
                return EncZonotopeStub {ez.mEncCenter, ez.mGenerators, ez.mP};
            }

            static inline ZonotopeStub stub(const Zonotope& z)
            {
                return ZonotopeStub {z.mCenter, z.mGenerators};
            }
        }
    }
}
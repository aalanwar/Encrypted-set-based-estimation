#include <EncEntity/EncEntity.h>

namespace EncEst {

EncEntity::EncEntity(const matrix_t<double>& HH, const matrix_t<double>& RR, Paillier* p):
                    mP(p),
                    mEncZ0(p->encryptVector_f(vector_t<double>::Zero(1,1)), matrix_t<double>::Zero(1,1), p)
{
  int q = RR.rows();
  int N = RR.cols();
  int n = HH.cols();
  //std::cout << HH.block(q, 0, q, n) << "\n" << RR.block(0, 1, q, 1) << "\n" << N << std::endl;
  for(int i = 0; i < N; i++) {
    matrix_t<double> t1 = HH.block(i*q, 0, q, n);
    vector_t<encnum> t2 = p->encryptVector_f(vector_t<double>::Zero(q, 1));
    vector_t<double> t3 = RR.block(0, i, q, 1);
    mEncStripList.push_back(EncStrip(t1, t2, t3, p));
  }
  mDist = std::uniform_real_distribution<double>(-1, 1);
}



/*void EncEntity::measure(const vector_t<double>& x) {
  int N = static_cast<int>(mEncStripList.size());
  int n = mEncStripList[0].mH.cols();
  int q = mEncStripList[0].mH.rows();
  matrix_t<double> P (2*q*N, n);
  vector_t<double> z (2*q*N, 1);
  P.setZero();
  z.setZero();
  for(int i = 0; i < N; i ++) {
    vector_t<double> noise (q, 1);
    for(int j = 0; j < q; j++)
      noise(j) = mEncStripList[i].mR(j)*mDist(mGen);
    mEncStripList[i].mEncY = mP.encryptVector_f(mEncStripList[i].mH*x + noise);
    P.block(2*i*q, 0, q, n) = mEncStripList[i].mH;
    P.block((2*i+1)*q, 0, q, n) = -mEncStripList[i].mH;
    z.block(2*i*q, 0, q, 1) = mEncStripList[i].mR + mEncStripList[i].mH*x + noise;
    z.block((2*i+1)*q, 0, q, 1) = mEncStripList[i].mR - (mEncStripList[i].mH*x + noise);
  }
  HPolytope<double> poly (P, z);
  auto zono = Converter<double>::toZonotope(poly);
  mEncZ0.mEncCenter = mP.encryptVector_f(zono.center());
  mEncZ0.mGenerators = zono.generators();
  mEncZ0.mDimension = zono.dimension();
  //std::cout << mEncZ0.mGenerators << std::endl;
}*/

}//end namespace

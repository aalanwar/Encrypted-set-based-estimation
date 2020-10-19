#include <EncStrip/EncStrip.h>
#include <stdexcept>
namespace EncEst {

EncStrip::EncStrip(const matrix_t<double>& H, const vector_t<encnum>& ency, const vector_t<double> r, Paillier* p):
 														mH(H),
														mEncY(ency),
														mR(r),
														mP(p)
														{}
EncStrip::EncStrip(const Strip& strip, Paillier* p):
                            mH(strip.mH),
                            mEncY(p->encryptVector_f(strip.mY)),
                            mR(strip.mR),
                            mP(p)
                            {}


}//end namespace

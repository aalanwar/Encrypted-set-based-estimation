#include <ConZonotope/ConZonotope.h>
#include <stdexcept>

namespace EncEst {

ConZonotope::ConZonotope(const vector_t<double>& center, const matrix_t<double>& generators, const matrix_t<double>& A, const vector_t<double>& b):
												mDimension(center.rows()),
												mCenter(center),
												mGenerators(generators),
												mA(A),
												mb(b)
												{}
ConZonotope ConZonotope::operator+(const ConZonotope& ez) {
	//cout << "mP =" << mP <<endl;
	//cout << "ez.mP =" << (ez.mP) <<endl; 
	if(mDimension != ez.mDimension)
		throw std::logic_error("Dimensions do not match!");
	matrix_t<double> G1 = mGenerators;
	matrix_t<double> G2 = ez.mGenerators;
	matrix_t<double> G (G1.rows(), G1.cols() + G2.cols());
	G << G1, G2;
	
	vector_t<double> newCenter =mCenter+ ez.mCenter;
	
	matrix_t<double> newA (mA.rows() + ez.mA.rows(), mA.cols() + ez.mA.cols());
	newA << mA, matrix_t<double>::Zero(mA.rows(),ez.mA.cols()),
			matrix_t<double>::Zero(ez.mA.rows(),mA.cols()), ez.mA;
	vector_t<double> newb (mb.rows() + ez.mb.rows(),mb.cols()) ;
	newb << mb, 
			ez.mb;
	return ConZonotope(newCenter, G,newA, newb);
}

ConZonotope ConZonotope::intersect(const ConZonotope& ey) {
	
	vector_t<double> newb (mb.rows() + ey.mb.rows() + mCenter.size(),mb.cols()) ;
	newb << mb, 
			   ey.mb,
			  ey.mCenter-mCenter;

	matrix_t<double> newA (mA.rows() + ey.mA.rows() + mGenerators.rows(), mA.cols() + ey.mA.cols());
	newA << mA, matrix_t<double>::Zero(mA.rows(),ey.mA.cols()),
			matrix_t<double>::Zero(ey.mA.rows(),mA.cols()), ey.mA,
			mGenerators, -ey.mGenerators;

	matrix_t<double> newG (mGenerators.rows(), newA.cols() );
	newG << mGenerators, matrix_t<double>::Zero(mGenerators.rows(),ey.mA.cols());		   
	return ConZonotope(mCenter, newG,newA, newb);
}

ConZonotope ConZonotope::intersect(const std::vector<Strip>& strips) {
	matrix_t<double> newA(mA.rows()+ strips.size(),mA.cols()+strips.size());
	newA.block(0,0,mA.rows(),newA.cols()) << mA,matrix_t<double>::Zero(mA.rows(),strips.size()); 
	
	matrix_t<double> newG (mGenerators.rows(),newA.cols());
	newG << mGenerators,matrix_t<double>::Zero(mGenerators.rows(),newA.cols() - mGenerators.cols()) ;

	vector_t<double> newb(mb.rows() + strips.size(), mb.cols());

	newb.block(0,0,mb.rows(),mb.cols()) = mb;

	//newA.block(mA.rows(),0,1,mA.cols()+strips.size()) << strips[i].mH * mGenerators , -strips[0].mR[0], matrix_t<double>::Zero(1,strips.size()-1);

	for (int i=0;i<strips.size();i++)
	{
		newA.block(mA.rows()+i,0,1,mA.cols()+strips.size()) << strips[i].mH * mGenerators , matrix_t<double>::Zero(1,i) , -strips[i].mR[0], matrix_t<double>::Zero(1,strips.size()-1-i);
		newb.block(i+1,0,1,1) =strips[i].mY - (strips[i].mH * mCenter);
	}
	//cout << newA << endl;
	/*matrix_t<double> tmpG2 (rstrips.size() , resG2.cols() + 1);
	for(auto i = 1; i < strips.size(); i++) {		
		matrix_t<double> tmpG2 (newA.rows()+1 , resG2.cols() + 1);
		tmpG2 << resG2, lambdas[i] * strips[i].mR;
		resG2 = tmpG2;
	}
	*/
	return ConZonotope(mCenter, newG,newA, newb );

}

}

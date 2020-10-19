#include "Zonotope.h"
#include <stdexcept>

namespace EncEst {

Zonotope::Zonotope(const vector_t<double>& center, const matrix_t<double>& generators):
												mDimension(center.rows()),
												mCenter(center),
												mGenerators(generators)
												{}



Zonotope Zonotope::intersect(const std::vector<Strip>& strips) {
	int n = mCenter.rows();
	int q = strips[0].mH.rows();
	int m = mGenerators.cols();
	std::vector<matrix_t<double>> lambdas= computeLambdas(strips);
	vector_t<double> resc = mCenter + lambdas[0]*(strips[0].mY - (strips[0].mH * mCenter));
	matrix_t<double> resG1 = mGenerators - (lambdas[0]*strips[0].mH * mGenerators);

	matrix_t<double> resG2 =  strips[0].mR[0] * lambdas[0] ;
	//matrix_t<double> resG2 =  55* lambdas[0] ;
	for(auto i = 1; i < strips.size(); i++) 
	{
		resc = resc  + lambdas[i]*(strips[i].mY - (strips[i].mH * mCenter)); 
		resG1 = resG1 - lambdas[i] * strips[i].mH * mGenerators;
		matrix_t<double> tmpG2 (resG2.rows(), resG2.cols() + 1);
		tmpG2 << resG2, lambdas[i] * strips[i].mR;
		resG2 = tmpG2;
	}

	matrix_t<double> resG (resG1.rows(), resG1.cols() + resG2.cols());
	
	resG << resG1, resG2;
	
	return Zonotope(resc, resG);
}

std::vector<matrix_t<double>> Zonotope::computeLambdas(const std::vector<Strip>& strips) {
	//
	//matrix_t<double> h_combined (strips[0].mH.rows() *strips.size() )( strips[0].mH.size.cols() );
	matrix_t<double> h_combined = strips[0].mH;

	for(int i=1; i< strips.size(); i++)
	{
		matrix_t<double> tmph (h_combined.rows()+ strips[i].mH.rows(), h_combined.cols() );
		tmph << h_combined, 
			    strips[i].mH;
		h_combined = tmph;				  
	}

	matrix_t<double> gamma = matrix_t<double>::Identity(strips.size(),strips.size());
	matrix_t<double> num = mGenerators * mGenerators.transpose() * h_combined.transpose();
	matrix_t<double> den = h_combined * mGenerators * mGenerators.transpose() * h_combined.transpose();
	for(int i=0; i< strips.size(); i++)
	{
		//Block of size (p,q), starting at (a,b) 
		//matrix.block(a,b,p,q);
		den = den + gamma.block(0,i,strips.size(),1) * strips[i].mR*strips[i].mR.transpose()*gamma.block(0,i,strips.size(),1).transpose();
	}

	//std::cout << h_combined<< endl;;
	//std::cout << h_combined[1];
	matrix_t<double> fulllambdas = num *(den.inverse());
	std::vector<matrix_t<double>> lambdas;
	//std::cout<<"lambda \n" << fulllambdas <<"\n";
	for(int i=0; i< strips.size(); i++)
	{
		lambdas.push_back( fulllambdas.block(0,i,fulllambdas.rows(),1) );
	}
	return lambdas;
}

}

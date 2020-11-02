#include <iostream>		// TODO: remove
#include <stdexcept>
#include <algorithm>

#include "Zonotope.h"

namespace EncEst 
{
	Zonotope::Zonotope(const vector_t<double>& center, const matrix_t<double>& generators)
		: mDimension(center.rows())
		, mCenter(center)
		, mGenerators(generators)
	{}

	Zonotope Zonotope::intersect(const std::vector<Strip>& strips) 
	{
		// TODO:SJ: add assertions here ...
		assert(std::none_of(strips.begin(), strips.end(), [this](const Strip& s) {return s.mH.cols() != mDimension;}));

		// TODO:SJ: This is not used ... (remove ???)
		int n = mCenter.rows();
		int q = strips[0].mH.rows();
		int m = mGenerators.cols();

		std::vector<matrix_t<double>> lambdas = computeLambdas(strips);

		// SJ: Formula (12) initial value
		vector_t<double> resc = mCenter + lambdas[0] * (strips[0].mY - (strips[0].mH * mCenter));

		// SJ: Formula (13) initial value (sum multiplied by generator)
		matrix_t<double> resG1 = mGenerators - (lambdas[0] * strips[0].mH * mGenerators);

		// SJ: Formula (13) initial value (lambda multplied by uncertainty)
		matrix_t<double> resG2 =  strips[0].mR[0] * lambdas[0];
		//matrix_t<double> resG2 =  55* lambdas[0] ;

		for(auto i = 1; i < strips.size(); i++) 
		{
			// SJ: Formula (12) expansion
			resc = resc  + lambdas[i] * (strips[i].mY - (strips[i].mH * mCenter));

			// SJ: Formula (13) expansion
			resG1 = resG1 - lambdas[i] * strips[i].mH * mGenerators;

			// SJ: Formla (13) append row to resG2
			matrix_t<double> tmpG2(resG2.rows(), resG2.cols() + 1);
			tmpG2 << resG2, lambdas[i] * strips[i].mR;
			resG2 = tmpG2;
		}

		// SJ: Concatenate resG1 and resG2.
		matrix_t<double> resG (resG1.rows(), resG1.cols() + resG2.cols());
		resG << resG1, resG2;
		
		return Zonotope(resc, resG);
	}

	std::vector<matrix_t<double>> Zonotope::computeLambdas(const std::vector<Strip>& strips) 
	{
		// TODO:SJ: Check if this method shouldn't be private/protected ... 

		//
		//matrix_t<double> h_combined (strips[0].mH.rows() *strips.size() )( strips[0].mH.size.cols() );
		
		// SJ: append all mH vertically of all strips
		matrix_t<double> h_combined = strips[0].mH;

		for(int i = 1; i < strips.size(); i++)
		{
			matrix_t<double> tmph (h_combined.rows() + strips[i].mH.rows(), h_combined.cols());
			tmph << h_combined, strips[i].mH;
			h_combined = tmph;				  
		}

		// TODO: remove
		std::cout << h_combined << std::endl;

		matrix_t<double> gamma = matrix_t<double>::Identity(strips.size(), strips.size());
		matrix_t<double> num = mGenerators * mGenerators.transpose() * h_combined.transpose();

		// SJ: Couldn't that be h_combined * num ????
		matrix_t<double> den = h_combined * mGenerators * mGenerators.transpose() * h_combined.transpose();

		for(int i = 0; i < strips.size(); i++)
		{
			//Block of size (p,q), starting at (a,b) 
			//matrix.block(a,b,p,q);
			// SJ: TODO, analyse why UnitTest failes here ...
			// SJ: gamma.block(...) are the e1,e2,... column-vectors of the identity matrix
			den = den + gamma.block(0, i, strips.size(), 1) * strips[i].mR * strips[i].mR.transpose() * gamma.block(0, i, strips.size(), 1).transpose();
		}

		matrix_t<double> fulllambdas = num * den.inverse();
		std::vector<matrix_t<double>> lambdas;

		for(int i = 0; i < strips.size(); i++)
		{
			// SJ: Add column-vectors to lambdas collection
			lambdas.push_back(fulllambdas.block(0, i, fulllambdas.rows(), 1));
		}

		return lambdas;
	}
}

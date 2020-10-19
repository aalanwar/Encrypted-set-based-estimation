#include <EncConZonotope/EncConZonotope.h>
#include <stdexcept>
//#include <glpk.h>
#include <eigen3/Eigen/Dense>
#include <eigen2mat/eigen2mat.h>
#include <vector>
using namespace std;
using namespace Eigen;
namespace EncEst {

EncConZonotope::EncConZonotope(const vector_t<encnum>& enccenter, const matrix_t<double>& generators, const matrix_t<double>& A,const vector_t<encnum>& encb,Paillier* p):
												mDimension(enccenter.rows()),
												mEncCenter(enccenter),
												mGenerators(generators),
												mA(A),
												mEncb(encb),
												mP(p)
												{}


EncConZonotope EncConZonotope::intersectMany(const std::vector<EncConZonotope>& conzonol) {
	EncConZonotope encconz_result = conzonol[0];
	for(auto i = 1; i < conzonol.size(); i++) 
	{
		encconz_result = encconz_result.intersect(conzonol[i]);
	}
	return encconz_result ;
}


EncConZonotope EncConZonotope::operator+(const EncConZonotope& ez) {
	//cout << "mP =" << mP <<endl;
	//cout << "ez.mP =" << (ez.mP) <<endl; 
	if(mP!=(ez.mP))
		throw std::logic_error("Public keys don't match!");
	else if(mDimension != ez.mDimension)
		throw std::logic_error("Dimensions do not match!");
	matrix_t<double> G1 = mGenerators;
	matrix_t<double> G2 = ez.mGenerators;
	matrix_t<double> G (G1.rows(), G1.cols() + G2.cols());
	G << G1, G2;
	
	vector_t<encnum> newEncCenter =mP->addVector_f(mEncCenter, ez.mEncCenter);
	
	matrix_t<double> newA (mA.rows() + ez.mA.rows(), mA.cols() + ez.mA.cols());
	newA << mA, matrix_t<double>::Zero(mA.rows(),ez.mA.cols()),
			matrix_t<double>::Zero(ez.mA.rows(),mA.cols()), ez.mA;
	vector_t<encnum> newEncb (mEncb.rows() + ez.mEncb.rows(),mEncb.cols()) ;
	newEncb << mEncb, 
			   ez.mEncb;
	return EncConZonotope(newEncCenter, G,newA, newEncb ,mP);
}

EncConZonotope EncConZonotope::intersect(const EncConZonotope& ey) {
	
	vector_t<encnum> newEncb (mEncb.rows() + ey.mEncb.rows() + mEncCenter.size(),mEncb.cols()) ;
	newEncb << mEncb, 
			   ey.mEncb,
			  mP->subVector_f(ey.mEncCenter,mEncCenter);

	matrix_t<double> newA (mA.rows() + ey.mA.rows() + mGenerators.rows(), mA.cols() + ey.mA.cols());
	newA << mA, matrix_t<double>::Zero(mA.rows(),ey.mA.cols()),
			matrix_t<double>::Zero(ey.mA.rows(),mA.cols()), ey.mA,
			mGenerators, -ey.mGenerators;

	matrix_t<double> newG (mGenerators.rows(), newA.cols() );
	newG << mGenerators, matrix_t<double>::Zero(mGenerators.rows(),ey.mA.cols());		   
	return EncConZonotope(mEncCenter, newG,newA, newEncb ,mP);
}

EncConZonotope EncConZonotope::intersect(const std::vector<EncStrip>& encstrips) {
	int n = mEncCenter.rows();
	int q = encstrips[0].mH.rows();
	int m = mGenerators.cols();
	Paillier* p = mP;
	matrix_t<double> newA(mA.rows()+ encstrips.size(),mA.cols()+encstrips.size());
	newA.block(0,0,mA.rows(),newA.cols()) << mA,matrix_t<double>::Zero(mA.rows(),encstrips.size()); 
	
	matrix_t<double> newG (mGenerators.rows(),newA.cols());
	newG << mGenerators,matrix_t<double>::Zero(mGenerators.rows(),newA.cols() - mGenerators.cols()) ;

	vector_t<encnum> newEncb(mEncb.rows() + encstrips.size(), mEncb.cols());

	newEncb.block(0,0,mEncb.rows(),mEncb.cols()) = mEncb;

	//newA.block(mA.rows(),0,1,mA.cols()+encstrips.size()) << encstrips[i].mH * mGenerators , -encstrips[0].mR[0], matrix_t<double>::Zero(1,encstrips.size()-1);

	for (int i=0;i<encstrips.size();i++)
	{
		newA.block(mA.rows()+i,0,1,mA.cols()+encstrips.size()) << encstrips[i].mH * mGenerators , matrix_t<double>::Zero(1,i) , -encstrips[i].mR[0], matrix_t<double>::Zero(1,encstrips.size()-1-i);
		newEncb.block(i+1,0,1,1) =mP->subVector_f(encstrips[i].mEncY,mP->multConstMatrixVector_f(encstrips[i].mH, mEncCenter));
	}
	//cout << newA << endl;
	/*matrix_t<double> tmpG2 (rencstrips.size() , resG2.cols() + 1);
	for(auto i = 1; i < encstrips.size(); i++) {		
		matrix_t<double> tmpG2 (newA.rows()+1 , resG2.cols() + 1);
		tmpG2 << resG2, lambdas[i] * encstrips[i].mR;
		resG2 = tmpG2;
	}
	*/
	return EncConZonotope(mEncCenter, newG,newA, newEncb ,mP);

}

/*
EncZonotope EncZonotope::intersectMany(const std::vector<EncZonotope>& zonol) {
	//if(weights.rows() != static_cast<int>(zonol.size()))
	//	throw std::logic_error("#weights and #EncZonotopes has to match!");
	int xx =1;
	vector_t<double> weights = computeWeights(zonol);
	//vector_t<double> weights = vector_t<double>::Ones(zonol.size(), 1);
	Paillier& p = zonol[0].mP;
	double sumW = weights.sum();
	auto res_c = p.multConstVector_f(weights(0)/sumW, zonol[0].mEncCenter);
	matrix_t<double> res_G = weights(0)/sumW * zonol[0].mGenerators;
	for(auto i = 1; i < zonol.size(); i++) {
		if(&p != &(zonol[i].mP))
			throw std::logic_error("Public keys don't match!");

		res_c = p.addVector_f(res_c, p.multConstVector_f(weights(i)/sumW, zonol[i].mEncCenter));
		auto Gtmp = zonol[i].mGenerators;
		matrix_t<double> G (res_G.rows(), res_G.cols() + Gtmp.cols());
		G << res_G, weights(i)/sumW * Gtmp;
		res_G = G;
	}
	return EncZonotope(res_c, res_G, p);
}



vector_t<double> EncZonotope::computeWeights(const std::vector<EncZonotope>& enczonos) {
	double sumofw =1;
	int e = static_cast<int>(enczonos.size());
	int n = enczonos[0].mGenerators.rows();
	matrix_t<double> A (e, e);
	double invtVecSum = 0;
	vector<double> tVec(e);
	vector_t<double> w(e);
	//vector_t<double> Gitmp;
	auto Gitmp = enczonos[0].mGenerators; //just for init
	for(int i = 0; i < e; i++) {
		Gitmp = enczonos[i].mGenerators;
		tVec[i] = (Gitmp*Gitmp.transpose()).trace();
		invtVecSum = invtVecSum + 1/tVec[i];
	}
	for(int i = 0; i < e; i++) {
		w[i] = sumofw/(tVec[i]*invtVecSum);
	}
	return w;
}

EncZonotope EncZonotope::intersect(const std::vector<EncStrip>& encstrips) {
	int n = mEncCenter.rows();
	int q = encstrips[0].mH.rows();
	int m = mGenerators.cols();
	std::vector<matrix_t<double>> lambdas= computeLambdas(encstrips);
	Paillier* p = mP;
	auto summand1 =mP->multConstMatrixVector_f(lambdas[0],mP->addVector_f(encstrips[0].mEncY,mP->multConstMatrixVector_f(-encstrips[0].mH, mEncCenter)));
	auto resc =mP->addVector_f(mEncCenter, summand1);
	matrix_t<double> resG1 = mGenerators - (lambdas[0]*encstrips[0].mH * mGenerators);

	cout << "\n lambdas0 = " <<lambdas[0] <<endl;
	cout << "\n lambdas1 = " <<lambdas[1] <<endl;
	cout << "\n mR = " << encstrips[0].mR <<endl;
	
	matrix_t<double> resG2 =  encstrips[0].mR[0] * lambdas[0] ;
	//matrix_t<double> resG2 =  55* lambdas[0] ;
	for(auto i = 1; i < encstrips.size(); i++) {
		if(&p != &(encstrips[i].mP))
			throw std::logic_error("Public keys don't match!");
		auto summandtmp =mP->multConstMatrixVector_f(lambdas[i],mP->addVector_f(encstrips[i].mEncY,mP->multConstMatrixVector_f(-encstrips[i].mH, mEncCenter)));
		resc =mP->addVector_f(resc, summandtmp);
		resG1 = resG1 - lambdas[i] * encstrips[i].mH * mGenerators;
		matrix_t<double> tmpG2 (resG2.rows(), resG2.cols() + 1);
		tmpG2 << resG2, lambdas[i] * encstrips[i].mR;
		resG2 = tmpG2;
	}
	matrix_t<double> resG (resG1.rows(), resG1.cols() + resG2.cols());
	resG << resG1, resG2;
	return EncZonotope(resc, resG, mP);
	//Zonotope<double> z (vector_t<double>::Zero(2,1), matrix_t<double>::Identity(2,2));
	//return EncZonotope(z, mP);
}

std::vector<matrix_t<double>> EncZonotope::computeLambdas(const std::vector<EncStrip>& encstrips) {
	//
	//matrix_t<double> h_combined (encstrips[0].mH.rows() *encstrips.size() )( encstrips[0].mH.size.cols() );
	matrix_t<double> h_combined = encstrips[0].mH;

	for(int i=1; i< encstrips.size(); i++)
	{
		matrix_t<double> tmph (h_combined.rows()+ encstrips[i].mH.rows(), h_combined.cols() );
		tmph << h_combined, 
			    encstrips[i].mH;
		h_combined = tmph;				  
	}
	matrix_t<double> gamma = matrix_t<double>::Identity(encstrips.size(),encstrips.size());
	matrix_t<double> num = mGenerators * mGenerators.transpose() * h_combined.transpose();
	matrix_t<double> den = h_combined * mGenerators * mGenerators.transpose() * h_combined.transpose();
	for(int i=0; i< encstrips.size(); i++)
	{
		//Block of size (p,q), starting at (a,b) 
		//matrix.block(a,b,p,q);
		den = den + gamma.block(0,i,encstrips.size(),1) * encstrips[i].mR*encstrips[i].mR.transpose()*gamma.block(0,i,encstrips.size(),1).transpose();
	}

	//std::cout << h_combined<< endl;;
	//std::cout << h_combined[1];
	matrix_t<double> fulllambdas = num *(den.inverse());
	std::vector<matrix_t<double>> lambdas;
	for(int i=0; i< encstrips.size(); i++)
	{
		lambdas.push_back( fulllambdas.block(0,i,encstrips.size(),1) );
	}
	return lambdas;
}
std::vector<matrix_t<double>> EncZonotope::computeLambdas(const std::vector<EncStrip>& encstrips) {
	int n = mGenerators.rows();
	int q = encstrips[0].mH.rows();
	int N = static_cast<int>(encstrips.size());
	auto GGt = mGenerators*mGenerators.transpose();
	matrix_t<double> A (N*n*q, N*n*q); //N strips, lambda_i collapsed n*q elements
	A.setZero();
	vector_t<double> b (N*n*q, 1);
	b.setZero();

	for(int i = 0; i < N; i++) {
		auto es = encstrips[i];
		auto C = GGt * es.mH.transpose();
		Matrix<double,Dynamic,Dynamic,RowMajor> Ctmp (C);
		Map<vector_t<double>> ci (Ctmp.data(), Ctmp.size());
		matrix_t<double> Hi (n*q, N*n*q);
		Hi.setZero();
		for(int j = 0; j < N; j++) {
			matrix_t<double> Htmp = encstrips[j].mH * GGt * es.mH.transpose();
			for(int k = 1; k <= n; k++) {
				Hi.block((k-1)*q, (k-1)*q + j*n*q, q, q) = Htmp.transpose();
			}
		}
		matrix_t<double> Ri (n*q, N*n*q);
		Ri.setZero();
		auto Rbt = es.mR*es.mR.transpose();
		for(int j = 1; j <= n; j++) {
			Ri.block((j-1)*q, (j-1)*q + i*n*q, q, q) = Rbt;
		}
		auto Ai = Ri + Hi;
		auto bi = ci;
		//Add rows
		A.block(i*n*q, 0, Ai.rows(), Ai.cols()) = Ai;
		b.block(i*n*q, 0, bi.rows(), 1) = bi;
	}
	matrix_t<double> AA (2*A.rows(),A.cols());
	AA << A,
			 -A;
	vector_t<double> bb (2*b.rows(),1);
	bb << b,
	      b;
	Optimizer<double> optimizer (AA, bb);
	EvaluationResult<double> sol = optimizer.getInternalPoint();
	vector_t<double> l = sol.optimumValue;
	std::vector<matrix_t<double>> lambdas;
	for(int i = 0; i < N; i++) {
		matrix_t<double> lambda = l.block(i*n*q, 0, n*q, 1);
		lambdas.push_back(Map<matrix_t<double>>(lambda.data(), n, q));
	}
	return lambdas;
}


	/*
	int e = static_cast<int>(enczonos.size());
	int n = enczonos[0].mGenerators.rows();
	matrix_t<double> A (e, e);
	for(int i = 0; i < e; i++) {
		for(int j = i; j < e; j++) {
			//pad smaller generator matrix to be able to multiply
			auto Gitmp = enczonos[i].mGenerators;
			auto Gjtmp = enczonos[j].mGenerators;
			int diff = Gitmp.cols() - Gjtmp.cols();
			if(diff > 0) {
				matrix_t<double> Gtmp (n, Gitmp.cols());
				Gtmp << Gjtmp, matrix_t<double>::Zero(n, diff);
				Gjtmp = Gtmp;
			}
			else if(diff < 0) {
				matrix_t<double> Gtmp (n, Gjtmp.cols());
				Gtmp << Gitmp, matrix_t<double>::Zero(n, -diff);
				Gitmp = Gtmp;
			}
			double tmp = (Gitmp*Gjtmp.transpose()).trace();
			A(i, j) = tmp;
			A(j, i) = tmp;
		}
	}
	matrix_t<double> Ainv = A.inverse();
	vector_t<double> one = vector_t<double>::Ones(e, 1);
	return Ainv*one / (one.transpose()*Ainv*one);
	


//global overloaded operators;
EncZonotope operator*(const matrix_t<double>& mat, EncZonotope &ez) {
	if(ez.mDimension != mat.cols())
		throw std::logic_error("Dimensions do not match!");
	return EncZonotope(ez.mP.multConstMatrixVector_f(mat, ez.mEncCenter),mat*ez.mGenerators, ez.mP);
}
*/
}

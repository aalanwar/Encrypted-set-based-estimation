/*
 * @file    example_zonotope.cpp
 * @author  Leonardo Winter Pereira
 * @since   2015-09-09
 * @version 2015-09-09
 */

#include <config.h>
#include <util/plotting/Plotter.h>
#include <representations/GeometricObject.h>
#include <sys/time.h>
#include <paillier.hh>
#include <gm.hh>
#include <NTL/ZZ.h>
#include <gmpxx.h>
#include <math/util_gmp_rand.h>
#include <math.h>
#include <ctime>
#include <assert.h>
#include <representations/EncZonotope/EncZonotope.h>
#include <eigen2mat/eigen2mat.h>
#include <map>
#include <random>
#include <representations/EncEntity/EncEntity.h>


//typedef int Number;
using namespace NTL;
using namespace hypro;
using namespace std;


int main()
{
	SetSeed(to_ZZ(time(NULL)));
	gmp_randstate_t randstate;
	gmp_randinit_default(randstate);
	gmp_randseed_ui(randstate,time(NULL));

	auto sk = Paillier_priv::keygen(randstate,1000,256); //600 ,256
	Paillier_priv pp(sk,randstate);

	auto pk = pp.pubkey();
	mpz_class nn = pk[0];
	Paillier p(pk,randstate);
	//////////
	int N = 100;
	matrix_t<double> F (2, 2);
	F << 0.992, -0.1247,
			 0.1247, 0.992;
	matrix_t<double> HH (4, 2);
	HH << 1, 0,
				0, 1,
				1, 0,
				0, 1;
	vector_t<double> rotc (2, 1);
	rotc << 5,
	        5;

	matrix_t<double> RR1 (1, 4);
	RR1 << 0.01, 0.05, 0.05, 0.02;
	matrix_t<double> RR2 (1, 3);
	RR2 << 0.1, 0.15, 0.1;
	matrix_t<double> RR3 (1, 2);
	RR3 << 0.02, 0.005;
	matrix_t<double> RR4 (1, 4);
	RR4 << 0.005, 0.2, 0.01, 0.05;
	matrix_t<double> RR5 (1, 2);
	RR5 << 0.01, 0.5;
	EncEntity encent1 (HH, RR1, p);
	EncEntity encent2 (HH.block(0,0,3,2), RR2, p);
	EncEntity encent3 (HH.block(0,0,2,2), RR3, p);
	EncEntity encent4 (HH, RR4, p);
	EncEntity encent5 (HH.block(0,0,2,2), RR5, p);
	vector<EncEntity> ents = {encent1, encent2, encent3, encent4, encent5};
	int e = static_cast<int>(ents.size());
	//since we are using hypro::converter, the resulting zonotope of halfspaces has m=n generators
	int n = HH.cols();
	int m = n;
	//maximum number of strips per entity
	int s = 4;
	matrix_t<double> resG ((e+1)*n*N, (m + s)*e);
	matrix_t<double> resc (n, (e+1)*N);
	resG.setZero();
	resc.setZero();
	matrix_t<double> xx (n,N);
	vector_t<double> x = rotc;
	for(int i = 0; i < N; i++) {
		vector<EncZonotope> enczonos;
		for(int j = 0; j < e; j++) {
			ents[j].measure(x);
			auto eztmp = ents[j].mEncZ0.intersect(ents[j].mEncStripList);
			resG.block(i*(e+1)*n + (j+1)*n, 0, n, eztmp.mGenerators.cols()) = eztmp.mGenerators;
			resc.block(0, i*(e+1)+j+1, n, 1) = pp.decryptVector_f(eztmp.mEncCenter);
			enczonos.push_back(eztmp);
		}
		auto weights = EncZonotope::computeWeights(enczonos);
		auto encres = EncZonotope::intersectMany(enczonos, weights);
		resG.block(i*(e+1)*n, 0, n, encres.mGenerators.cols()) = encres.mGenerators;
		resc.block(0, (e+1)*i, n, 1) = pp.decryptVector_f(encres.mEncCenter);

		xx.block(0, i, n, 1) = x;
		//Update position
		x = F*x;
	}
	matrix_t<double> etmp (1,1);
	etmp << e;
	map<const char*, const matrix_t<double>> varMap {
		{ "GG", resG },
    { "cc", resc },
		{ "xx", xx },
		{ "e", etmp }
	};
	Eigen2Mat::writeToFile(varMap, "amr_2.mat");
	cout << "\nworks!\n";
	/*
	matrix_t<double> G1 (3, 5);
	G1 << 1, 0, 1,-1, 0,
       -1,-1, 0,-2, 0,
        1,-1, 0, 0, 0;
	matrix_t<double> G2 (3, 5);
	G2 << 1,-1, 1, 0, 0,
       -1, 0,-1, 0, 0,
       -1, 3, 1, 0, 0;
	matrix_t<double> G3 (3, 5);
	G3 << 0, 1,-2, 0, 0,
        1, 2, 0, 0,-1,
       -2, 1, 0, 0,-1;
	EncZonotope ez1 (p.encryptVector_f(vector_t<double>::Constant(3, 2)), G1, p);
	EncZonotope ez2 (p.encryptVector_f(vector_t<double>::Constant(3, 3)), G2, p);
	EncZonotope ez3 (p.encryptVector_f(vector_t<double>::Constant(3, 4)), G3, p);
	vector<EncZonotope> ezs = {ez1, ez2, ez3};
	auto weights = EncZonotope::computeWeights(ezs);
	auto er = EncZonotope::intersectMany(ezs, weights);
	cout << "\n\nGenerators:\n" << er.mGenerators << "\nCenter:\n" << pp.decryptVector_f(er.mEncCenter) << "\nWeights:\n"<< weights << "\n\n";
	*/
  return 0;
}

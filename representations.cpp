#include <sys/time.h>
#include "crypto/paillier.hh"
#include <gm.hh>
#include <NTL/ZZ.h>
#include <gmpxx.h>
#include <math/util_gmp_rand.h>
#include <math.h>
#include <ctime>
#include <assert.h>
#include <EncZonotope/EncZonotope.h>
#include <Zonotope/Zonotope.h>
#include <EncConZonotope/EncConZonotope.h>
#include <ConZonotope/ConZonotope.h>
#include <eigen2mat/eigen2mat.h>
#include <CSVReader.h>
#include <map>
#include <random>
//#include <representations/EncEntity/EncEntity.h>
#include <iostream>

using namespace NTL;
using namespace std;
using namespace EncEst;
using namespace Eigen;

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
	/*
	matrix_t<double> G1 (2, 2);
	G1 << 1, 3,
		  3, 5;
	vector_t<double> c1 (2,1);
	c1 << -1, 1;
	matrix_t<double> G2 (2, 3);
	G2 << 1, 2, 1,
		 -3, 4, 1;
	vector_t<double> c2 (2,1);
	c2 << -0.5, 1.5;
	Zonotope zono1 (c1,G1);
	Zonotope zono2 (c2,G2);
	EncZonotope encz1 = p.encrypt(zono1);
	EncZonotope encz2 = p.encrypt(zono2);
	auto res_intersect = encz1.intersect(encz2);
	auto res = pp.decrypt(res_intersect);
	std::cout << res.mCenter << "\n\n" << res.mGenerators << "\n\n";
	*/
	//////////////////////////////////
	/*
	matrix_t<double> h1(1,2);
	h1 << 1,0;
	vector_t<double> y1(1) ;
	y1 << -2;
    vector_t<double> R1(1) ;
	R1 << 5;
	Strip s1 (h1,y1,R1);
	EncStrip encs1 = p.encrypt(s1);

	matrix_t<double> h2(1,2);
	h2 << 0,1;
	vector_t<double> y2(1) ;
	y2 << 2;
    vector_t<double> R2(1) ;
	R2 << 3;
	Strip s2 (h2,y2,R2);
	EncStrip encs2 = p.encrypt(s2); 

	auto ss = pp.decrypt(encs1);
	std::cout << ss.mH << "\n\n" << ss.mY << "\n\n" << ss.mR << "\n\n";


	std::vector<EncStrip> encstrips={encs1,encs2};

	auto zonoStrEnc = encz1.intersect(encstrips);
	std::cout << "\n Zonotope Strips Intersections \n";
	auto zonoStr = pp.decrypt(zonoStrEnc);
	std::cout << zonoStr.mCenter << "\n\n" << zonoStr.mGenerators << "\n\n";
	*/
    ///////////////////////////////////////
	vector_t<double> c1 (2,1);
	c1 << 0,
		  0;
	matrix_t<double> G1 (2, 3);
	G1 << 3, 0,1,
		  0, 2,1;

	matrix_t<double> A1 (1, 3);
	A1 << 1, 0, 1;
	vector_t<double> b1 (1,1);
	b1 << 1;
	ConZonotope conzono1 (c1,G1,A1,b1);
	EncConZonotope encconz1 = p.encrypt(conzono1);

	vector_t<double> c2 (2,1);
	c2 << 0,
		  0;
	matrix_t<double> G2 (2, 3);
	G2 << 1.5, -1.5,0.5,
		  1, 0.5,-1;

	matrix_t<double> A2 (1, 3);
	A2 << 1, 1, 1;
	vector_t<double> b2 (1,1);
	b2 << 1;
	ConZonotope conzono2 (c2,G2,A2,b2);




	EncConZonotope encconz2 = p.encrypt(conzono2);
	auto res_intersect = encconz1.intersect(encconz2);
	auto res = pp.decrypt(res_intersect);
	std::cout << res.mCenter << "\n\n" << res.mGenerators << "\n\n"<< res.mA << "\n\n" << res.mb << endl;
	std::cout << "\n ********************************** \n" << endl;

	matrix_t<double> h1(1,2);
	h1 << 1,0;
	vector_t<double> y1(1) ;
	y1 << -2;
    vector_t<double> R1(1) ;
	R1 << 5;
	Strip s1 (h1,y1,R1);
	EncStrip encs1 = p.encrypt(s1);

	matrix_t<double> h2(1,2);
	h2 << 0,1;
	vector_t<double> y2(1) ;
	y2 << 2;
    vector_t<double> R2(1) ;
	R2 << 3;
	Strip s2 (h2,y2,R2);
	EncStrip encs2 = p.encrypt(s2); 

	std::vector<EncStrip> encstrips={encs1,encs2};

	auto conzonoStrEnc = encconz1.intersect(encstrips);
	auto rescon = pp.decrypt(conzonoStrEnc);
	std::cout << rescon.mCenter << "\n\n" << rescon.mGenerators << "\n\n"<< rescon.mA << "\n\n" << rescon.mb << endl;
	//******************************************//
		// Creating an object of CSVWriter
	CSVReader reader("zvector.csv");
 
	// Get the data from CSV File
	std::vector<std::vector<std::string> > dataList = reader.getData();
 
	// Print the content of row by row on screen
	for(std::vector<std::string> vec : dataList)
	{
		//for(std::string data : vec)
		//{
		//	std::cout<<std::stof(data) << " , ";
		//}
		std::cout << std::stof(vec[0]);
		std::cout<<std::endl;
	}
	
	/*
	map<const char*, const matrix_t<double>> varMap {
		{ "Gres", res.mGenerators },
		{ "cres", res.mCenter }
	};
	
	Eigen2Mat::writeToFile(varMap, "amr_test.mat");
	*/
  return 0;
}
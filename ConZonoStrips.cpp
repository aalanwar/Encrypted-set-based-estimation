#include <sys/time.h>
#include "crypto/paillier.hh"
#include <gm.hh>
#include <NTL/ZZ.h>
#include <gmpxx.h>
#include <math/util_gmp_rand.h>
#include <algorithm>    // std::rotate
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


	//auto sk = Paillier_priv::keygen(randstate,1000,256); //600 ,256
	auto sk = Paillier_priv::keygen(randstate,1000,256);
	Paillier_priv pp(sk,randstate);

	auto pk = pp.pubkey();
	mpz_class nn = pk[0];
	Paillier p(pk,randstate);
	//////////
	
	vector_t<double> c1 (3,1);
	c1 << 0,
		  0,
		  0;
	matrix_t<double> G1 (3, 3);
	G1 << 7.5, 0,0,
		  0, 7.5,0,
		  0, 0  ,7.5;

	matrix_t<double> A1 (1, 3);
	A1 << 0, 0, 0;
	vector_t<double> b1 (1,1);
	b1 << 0;
	ConZonotope conzono (c1,G1,A1,b1);
	EncConZonotope encconz = p.encrypt(conzono);
	EncConZonotope encconz_init = encconz;
	
	vector_t<double> cQ (3,1);
	cQ << 0,
		  0,
		  0;
	matrix_t<double> GQ (3, 3);
	GQ << 0.1, 0,0,
		  0, 0.1,0,
		  0, 0  ,0.1;

	matrix_t<double> AQ (1, 3);
	AQ << 0, 0, 0;
	vector_t<double> bQ (1,1);
	bQ << 0;
	ConZonotope Qconzono (cQ,GQ,AQ,bQ);
	EncConZonotope Qencconz = p.encrypt(Qconzono);
	
	//std::cout << "\n Zonotope Strips Intersections \n";

	
    
	
	// Creating an object of CSVWriter
	//CSVReader reader("zvector_zcontstrips.csv",",");
	CSVReader reader("zvector.csv",",");
 
	// Get the data from CSV File
	std::vector<std::vector<std::string> > dataList = reader.getData();
 
	// Print the content of row by row on screen
	std::vector<double> z(dataList.size(),1);
	int index=0;
	for(std::vector<std::string> vec : dataList)
	{
		//each line of csv file is a vector (vec) of strings
		z[index++] = std::stof(vec[0]);
		
	}

	
	
	matrix_t<double> h(1,3);
	vector_t<double> y(1) ;
    vector_t<double> R(1) ;

	//for (int i=0; i <dataList.size();i=i+8)
	int fileindex = 1;
	int steps =10; //dataList.size();
	//std::vector<map> vecofmaps (steps);
	int matindex = 0;
	int sindex =0;
	map<const char*, const matrix_t<double>> varMap ;
	std::vector<string> s(dataList.size());
	clock_t start_time, end_time;
	std::vector<double> agg_time (dataList.size());
	int agg_time_index =0;
	std::vector<double> qry_time (dataList.size());
	int qry_time_index =0;
	for (int i=0; i < dataList.size()/3; i=i+3)
	//for (int i=0; i < 400; i=i+3)
	{
		std::vector<double> rotatingh = {1,0,0} ;
		std::vector<EncStrip> encstripsvec;
		encstripsvec.clear();
		for (int j=0; j<3; j++)
		{
			std::rotate(rotatingh.rbegin(),rotatingh.rbegin()+1,rotatingh.rend());
			h << rotatingh[0],rotatingh[1],rotatingh[2];
			//std::cout << h << endl;
			if(rotatingh[1]==1)
			{
				R << 0.8;
			}
			else
			{
				R << 0.2;
			}
			y << z[i+j];
			//cout<< y <<endl;
			Strip s (h,y,R);
			EncStrip encs  =  p.encrypt(s);
			encstripsvec.push_back(encs);
		}
		// to avoid overflow decrypt and encrypt
		//ConZonotope unencconz = pp.decrypt(encconz);
		//encconz =  p.encrypt(unencconz);
		start_time = clock();
		EncConZonotope res_intersect = encconz_init.intersect(encstripsvec);
		encconz = res_intersect+ Qencconz;
		end_time = clock();
		agg_time[agg_time_index++] = double(end_time - start_time)/ double(CLOCKS_PER_SEC);;

		start_time = clock();
		auto conzonoStr = pp.decrypt(encconz);
		end_time = clock();
		qry_time[qry_time_index++] = double(end_time - start_time)/ double(CLOCKS_PER_SEC);;

		//std::cout << "center =\n"<< conzonoStr.mCenter << "\n\n" << "Gen =\n"<< conzonoStr.mGenerators << "\n\n";
		//std::cout << "A =\n"<< conzonoStr.mA << "\n\n" << "b =\n"<< conzonoStr.mb << "\n\n";

		//std::cout<<i<<std::endl;

		s[sindex] = "c_"+ std::to_string(matindex);	
		varMap.insert(make_pair( s[sindex].c_str(), conzonoStr.mCenter ));

		s[sindex+1] = "G_"+ std::to_string(matindex);		
		varMap.insert(make_pair( s[sindex+1].c_str(), conzonoStr.mGenerators ));

		s[sindex+2] = "A_"+ std::to_string(matindex);		
		varMap.insert(make_pair( s[sindex+2].c_str(), conzonoStr.mA ));

		s[sindex+3] = "b_"+ std::to_string(matindex);		
		varMap.insert(make_pair( s[sindex+3].c_str(), conzonoStr.mb ));
		sindex = sindex+4;
		matindex  = matindex +1;

	}
	double avg = 0;
	for (int i=0; i<agg_time_index;i++ )
	{
		avg = avg + agg_time[i];
	}
	avg = avg/(agg_time_index-1);
	std::cout<<" agg time (sec)= "<< avg << endl;

	avg = 0;
	for (int i=0; i<qry_time_index;i++ )
	{
		avg = avg + qry_time[i];
	}
	avg = avg/(qry_time_index-1);
	std::cout<<" qry time (sec)= "<< avg << endl;

		std::string filename = "MATLAB/CMatFiles/cppCZonoStrips.mat";		
		Eigen2Mat::writeToFile(varMap, filename.c_str());
	/*
	map<const char*, const matrix_t<double>> varMap {
		{ "Gres", res.mGenerators },
		{ "cres", res.mCenter }
	};
	
	Eigen2Mat::writeToFile(varMap, "amr_test.mat");
	*/
  return 0;
}
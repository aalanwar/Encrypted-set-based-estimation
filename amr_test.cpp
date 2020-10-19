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
#include <Strip/Strip.h>
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
	ConZonotope conz_init = conzono;
	EncConZonotope encconz = p.encrypt(conzono);
	EncConZonotope encconz_init = encconz;
	

	
	ConZonotope conz_e1 = conzono;
	EncConZonotope encconz_e1 = encconz_init;

	ConZonotope conz_e2 = conzono;
	EncConZonotope encconz_e2 = encconz_init;

	ConZonotope conz_e3 = conzono;
	EncConZonotope encconz_e3 = encconz_init;
	//std::cout << "\n Zonotope Strips Intersections \n";


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

	
	
    
	
	// Creating an object of CSVWriter
	CSVReader reader("zvector_ent.csv",",");
 
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
	int zindex = 0;
	int matindex = 0;
	int sindex =0;
	map<const char*, const matrix_t<double>> varMap ;
	int steps = 4; //dataList.size()
	std::vector<string> s (dataList.size());
	clock_t start_time, end_time;
	std::vector<double> sensor_time (dataList.size());
	int sensor_time_index =0;
	std::vector<double> agg_time (dataList.size());
	int agg_time_index =0;
	std::vector<double> qry_time (dataList.size());
	int qry_time_index =0;
	for (int i=0; i < dataList.size()/(8*3); i++)
	//for (int i=0; i < 4; i++)
	{
		std::vector<Strip> stripsvec;
		// -----------------  e1  ------------------------
		std::vector<double> rotatingh = {1,0,0} ;
		stripsvec.clear();
		for (int j=0; j<3; j++)
		{
			std::rotate(rotatingh.rbegin(),rotatingh.rbegin()+1,rotatingh.rend());
			h << rotatingh[0],rotatingh[1],rotatingh[2];
			//std::cout << h << endl;
			if(rotatingh[1]==1)
			{
				R << 0.92;
			}
			else
			{
				R << 0.2;
			}
			y << z[zindex];
			zindex = zindex +1;
			//cout<< y <<endl;
			Strip s (h,y,R);
			
			stripsvec.push_back(s);
		}
		start_time = clock();
		conz_e1 = conz_e1.intersect(stripsvec);
		encconz_e1 =  p.encrypt(conz_e1);
		end_time = clock();
		sensor_time[sensor_time_index++] = double(end_time - start_time)/ double(CLOCKS_PER_SEC);;

		// -----------------  e2  ------------------------
		rotatingh = {1,0,0} ;
		stripsvec.clear();
		for (int j=0; j<3; j++)
		{
			std::rotate(rotatingh.rbegin(),rotatingh.rbegin()+1,rotatingh.rend());
			h << rotatingh[0],rotatingh[1],rotatingh[2];
			//std::cout << h << endl;
			if(rotatingh[1]==1)
			{
				R << 0.92;
			}
			else
			{
				R << 0.2;
			}
			y << z[zindex];
			zindex = zindex +1;
			//cout<< y <<endl;
			Strip s (h,y,R);
			
			stripsvec.push_back(s);
		}
		
		conz_e2 = conz_e2.intersect(stripsvec);
		encconz_e2 =  p.encrypt(conz_e2);
		// -----------------  e3  ------------------------
		rotatingh = {1,0,0} ;
		stripsvec.clear();
		for (int j=0; j<2; j++)
		{
			std::rotate(rotatingh.rbegin(),rotatingh.rbegin()+1,rotatingh.rend());
			h << rotatingh[0],rotatingh[1],rotatingh[2];
			//std::cout << h << endl;
			if(rotatingh[1]==1)
			{
				R << 0.92;
			}
			else
			{
				R << 0.2;
			}
			y << z[zindex];
			zindex = zindex +1;
			//cout<< y <<endl;
			Strip s (h,y,R);
			
			stripsvec.push_back(s);
		}
		
		conz_e3 = conz_e3.intersect(stripsvec);
		encconz_e3 =  p.encrypt(conz_e3);
		// ---------------- Aggregator ---------------------//
		start_time = clock();
		std::vector<EncConZonotope> conzonol = {encconz_e1, encconz_e2,encconz_e3};
		encconz = EncConZonotope::intersectMany(conzonol);

		// Time update
		encconz = encconz + Qencconz;
		end_time = clock();
		agg_time[agg_time_index++] = double(end_time - start_time)/ double(CLOCKS_PER_SEC);;


		auto unenczono = pp.decrypt(encconz);
		conz_e1 = conz_init;
		conz_e2 = conz_init;
		conz_e3 = conz_init;

		//std::cout << "center =\n"<< unenczono.mCenter << "\n\n" << "Gen =\n"<< unenczono.mGenerators << "\n\n";
		//std::cout << "A =\n"<< unenczono.mA << "\n\n" << "b =\n"<< unenczono.mb << "\n\n";
		//std::cout<<i<<std::endl;	
		s[sindex] = "c_"+ std::to_string(matindex);	
		varMap.insert(make_pair( s[sindex].c_str(), unenczono.mCenter ));

		s[sindex+1] = "G_"+ std::to_string(matindex);		
		varMap.insert(make_pair( s[sindex+1].c_str(), unenczono.mGenerators ));

		s[sindex+2] = "A_"+ std::to_string(matindex);		
		varMap.insert(make_pair( s[sindex+2].c_str(), unenczono.mA ));

		s[sindex+3] = "b_"+ std::to_string(matindex);		
		varMap.insert(make_pair( s[sindex+3].c_str(), unenczono.mb ));
		sindex = sindex+4;
		matindex  = matindex +1;

	}
		std::string filename = "MATLAB/CMatFiles/cppCZonoEntities.mat";		
		Eigen2Mat::writeToFile(varMap, filename.c_str());
		double avg = 0;
	for (int i=0; i<sensor_time_index;i++ )
	{
		avg = avg + sensor_time[i];
	}
	avg = avg/(sensor_time_index-1);
	std::cout<<" sensor time (sec)= "<< avg << endl;
	 
	avg = 0;
	for (int i=0; i<agg_time_index;i++ )
	{
		avg = avg + agg_time[i];
	}
	avg = avg/(agg_time_index-1);
	std::cout<<" agg time (sec)= "<< avg << endl;

	/*
	map<const char*, const matrix_t<double>> varMap {
		{ "Gres", res.mGenerators },
		{ "cres", res.mCenter }
	};
	
	Eigen2Mat::writeToFile(varMap, "amr_test.mat");
	*/
  return 0;
}
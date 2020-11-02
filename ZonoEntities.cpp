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
	
	matrix_t<double> G (3, 3);
	G << 7.5,0 ,0,
		  0 ,7.5,0,
		  0 ,0 ,7.5;

	vector_t<double> c (3,1);
	c << 0,0,0;

	Zonotope zono (c,G);
	EncZonotope enczono = p.encrypt(zono);
	
	Zonotope zono_e1 (c,G);
	EncZonotope enczono_e1 = p.encrypt(zono_e1);

	Zonotope zono_e2 (c,G);
	EncZonotope enczono_e2 = p.encrypt(zono_e2);

	Zonotope zono_e3 (c,G);
	EncZonotope enczono_e3 = p.encrypt(zono_e3);
	//std::cout << "\n Zonotope Strips Intersections \n";

	
	// Creating an object of CSVWriter
	CSVReader reader("zvector_ent.csv",",");
 
	// Get the data from CSV File
	std::vector<std::vector<std::string>> dataList = reader.getData();
 
	// Print the content of row by row on screen
	std::vector<double> z(dataList.size(), 1);
	int index=0;
	clock_t start_time, end_time;
	
	std::vector<double> sensor_time (dataList.size());
	int sensor_time_index = 0;
	
	std::vector<double> agg_time (dataList.size());
	int agg_time_index = 0;
	
	std::vector<double> qry_time (dataList.size());
	int qry_time_index = 0;

	for(std::vector<std::string> vec : dataList)
	{
		//each line of csv file is a vector (vec) of strings
		z[index++] = std::stof(vec[0]);
	}

	matrix_t<double> h(1,3);
	vector_t<double> y(1);
    vector_t<double> R(1);
	matrix_t<double> Q(3,3);
	Q << 0.1, 0 ,0,
		 0  ,0.1,0,
		 0  ,0  ,0.1;

	//for (int i=0; i <dataList.size();i=i+8)
	int fileindex = 1;
	int zindex = 0;
	//------------------ for converting to Matlab
	map<const char*, const matrix_t<double>> varMap ;
	int steps =10;//dataList.size();
	std::vector<string> s (dataList.size());
	int matindex = 0;
	int sindex =0;
	//-------------------
	for (int i = 0; i < dataList.size() / 3; i = i + 8)
	//for (int i=0; i < 5; i++)
	{
		std::vector<Strip> stripsvec;
		// -----------------  e1  ------------------------
		std::vector<double> rotatingh = {1,0,0};
		stripsvec.clear();

		for (int j = 0; j < 3; j++)
		{
			std::rotate(rotatingh.rbegin(), rotatingh.rbegin() + 1, rotatingh.rend());
			h << rotatingh[0], rotatingh[1], rotatingh[2];
			//std::cout << h << endl;
			if(rotatingh[1] == 1)
			{
				R << 0.78;
			}
			else
			{
				R << 0.15;
			}

			y << z[zindex];
			zindex = zindex +1;
			//cout<< y <<endl;
			Strip s (h,y,R);
			
			stripsvec.push_back(s);
		}
		// to avoid overflow decrypt and encrypt
		start_time = clock(); 
		zono_e1 = zono_e1.intersect(stripsvec);
		enczono_e1 =  p.encrypt(zono_e1);
		end_time = clock();
		sensor_time[sensor_time_index++] = double(end_time - start_time) / double(CLOCKS_PER_SEC);

		// -----------------  e2  ------------------------
		rotatingh = {1,0,0};
		stripsvec.clear();

		for (int j = 0; j < 3; j++)
		{
			std::rotate(rotatingh.rbegin(),rotatingh.rbegin()+1,rotatingh.rend());
			h << rotatingh[0],rotatingh[1],rotatingh[2];
			//std::cout << h << endl;
			if(rotatingh[1] == 1)
			{
				R << 0.78;
			}
			else
			{
				R << 0.15;
			}

			y << z[zindex];
			zindex = zindex +1;
			//cout<< y <<endl;
			Strip s (h,y,R);
			
			stripsvec.push_back(s);
		}

		// to avoid overflow decrypt and encrypt
		zono_e2 = zono_e2.intersect(stripsvec);
		enczono_e2 = p.encrypt(zono_e2);

		// -----------------  e3  ------------------------
		rotatingh = {1,0,0} ;
		stripsvec.clear();

		for (int j=0; j<2; j++)
		{
			std::rotate(rotatingh.rbegin(),rotatingh.rbegin()+1,rotatingh.rend());
			h << rotatingh[0],rotatingh[1],rotatingh[2];
			//std::cout << h << endl;
			if(rotatingh[1] == 1)
			{
				R << 0.78;
			}
			else
			{
				R << 0.15;
			}

			y << z[zindex];
			zindex = zindex +1;
			//cout<< y <<endl;
			Strip s (h,y,R);
			
			stripsvec.push_back(s);
		}

		// to avoid overflow decrypt and encrypt
		zono_e3 = zono_e3.intersect(stripsvec);
		enczono_e3 =  p.encrypt(zono_e3);
		
		// ---------------- Aggregator ---------------------//
		start_time = clock();
		std::vector<EncZonotope> zonol = {enczono_e1, enczono_e2,enczono_e3};
		enczono = enczono_e1.intersectMany(zonol);

		// Time update
		matrix_t<double> tempGen(enczono.mGenerators.rows(),enczono.mGenerators.cols()+Q.cols());
		tempGen.block(0,0,enczono.mGenerators.rows(),enczono.mGenerators.cols()) = enczono.mGenerators;
		tempGen.block(0,enczono.mGenerators.cols(),Q.rows(),Q.cols()) = Q;
		EncZonotope timeupdatezono (enczono.mEncCenter,tempGen,&p);
		enczono = timeupdatezono;
		enczono  = enczono.reduce(400);
		end_time = clock();
		agg_time[agg_time_index++] = double(end_time - start_time)/ double(CLOCKS_PER_SEC);;

		auto unenczono = pp.decrypt(enczono);
		enczono =  p.encrypt(unenczono);
		// for next iteration
		zono_e1 = unenczono;
		zono_e2 = unenczono;
		zono_e3 = unenczono;
		//////auto zonodebug = pp.decrypt(enczono);
		//std::cout << "center =\n"<< unenczono.mCenter << "\n\n" << "Gen =\n"<< unenczono.mGenerators << "\n\n";
		//std::cout<< i <<std::endl;
			
		
		s[sindex] = "c_"+ std::to_string(matindex);	
		varMap.insert(make_pair(s[sindex].c_str(), unenczono.mCenter));

		s[sindex+1] = "G_"+ std::to_string(matindex);		
		varMap.insert(make_pair(s[sindex+1].c_str(), unenczono.mGenerators));
		sindex = sindex + 2;
		matindex = matindex + 1;
	}

	double avg = 0;
	for (int i = 0; i < sensor_time_index; i++)
	{
		avg = avg + sensor_time[i];
	}

	avg = avg / (sensor_time_index - 1);
	std::cout<<" sensor time (sec)= "<< avg << endl;
	 
	avg = 0;
	for (int i = 0; i < agg_time_index; i++)
	{
		avg = avg + agg_time[i];
	}

	avg = avg / (agg_time_index - 1);
	std::cout<<" agg time (sec)= "<< avg << endl;

	std::string filename = "MATLAB/CMatFiles/cppZonoEntites.mat";		
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
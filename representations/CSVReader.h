#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
 
/*
 * A class to read data from a csv file.
 */
namespace EncEst {
class CSVReader
{
	std::string fileName;
	std::string delimeter;
 
public:
	CSVReader(std::string filename, std::string  delm) ;
	// Function to fetch data from a CSV File
	std::vector<std::vector<std::string> > getData();
};
 
}
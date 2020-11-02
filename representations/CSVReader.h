#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace EncEst 
{
	/// @brief This is some class description ....
	class CSVReader
	{
	private:
		std::string fileName;
		std::string delimeter;
	
	public:
		/// @brief Reads a CSV file and splits the lines by delm.
		/// @param filename the full path to the CSV file.
		/// @param delm the delimiter used to split each line in the file.
		CSVReader(std::string filename, std::string  delm);

		/// @brief Get the data of the file split by the specified delimiter.
		/// @return vector of vector of std::strings which are the lines separated by the sepcified delimiter.
		std::vector<std::vector<std::string>> getData();

	};
}

#include <CSVReader.h>

namespace EncEst 
{
	CSVReader::CSVReader(std::string filename, std::string delm = ",") 
		: fileName(filename)
		, delimeter(delm)
	{ }

	std::vector<std::vector<std::string>> CSVReader::getData()
	{
		std::ifstream file(fileName);
	
		std::vector<std::vector<std::string>> dataList;
	
		std::string line = "";

		while (getline(file, line))
		{
			std::vector<std::string> vec;
			boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
			dataList.push_back(vec);
		}

		file.close();
	
		return dataList;
	}

}
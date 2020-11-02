#pragma once

#include <mat.h>
#include <eigen3/Eigen/Dense>
#include <map>

class Eigen2Mat 
{
public:
  static int writeToFile(const std::map<const char*, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>& varMap, const char* file);
};

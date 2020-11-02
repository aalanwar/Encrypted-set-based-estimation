/*
 * matcreat.cpp - MAT-file creation program
 *
 * See the MATLAB External Interfaces/API Guide for compiling information.
 *
 * Calling syntax:
 *
 *   matcreat
 *
 * Create a MAT-file which can be loaded into MATLAB.
 *
 * This program demonstrates the use of the following functions:
 *
 *  matClose
 *  matGetVariable
 *  matOpen
 *  matPutVariable
 *  matPutVariableAsGlobal
 *
 * Copyright 1984-2007 The MathWorks, Inc. (modified)
 */

#include <stdlib.h>
#include <iostream>
#include "eigen2mat.h"

using namespace Eigen;

int Eigen2Mat::writeToFile(const std::map<const char*, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>& varMap, const char* file) 
{  
  MATFile *pmat;
  mxArray *pa;

  pmat = matOpen(file, "w");

  if (pmat == NULL) 
  {
    printf("Error creating file %s\n", file);
    printf("(Do you have write permission in this directory?)\n");
    return(EXIT_FAILURE);
  }

  for(const auto& [key, val] : varMap) 
  {
    auto M = val;
    const char* name = key;
    const double* data = M.data();
    int n = M.rows();
    int m = M.cols();

    pa = mxCreateDoubleMatrix(n, m, mxREAL);
    
    if (pa == NULL) 
    {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }

    memcpy((void *)(mxGetPr(pa)), (void *)data, sizeof(double) * m * n);

    int status = matPutVariable(pmat, name, pa);
    if (status != 0) 
    {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

    mxDestroyArray(pa);
  }

  if (matClose(pmat) != 0) 
  {
    printf("Error closing file %s\n",file);
    return(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
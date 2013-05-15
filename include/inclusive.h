#ifndef definitions_include_declared
#define definitions_include_declared
// release version 3.1.5 : see define VERSION

#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <cstring>
#include <cstdlib> // 
#include <algorithm> // 

using namespace std;

#define READ 0
#define WRITE 1

#define SINGLE 0
#define BOTH 1

// (3.1.2) 3 bugs fixed in Utilities.cpp 
// (3.1.2) 1 bug fixed in RandomNumber.cpp
// (3.1.2) 1 bug fixed in MotifSamplerRun.cpp (LL)
// (3.1.2) 1 bug fixed error input -x parameter
// (3.1.4) 1 bug fixed in MotifSamplerRun.cpp (mask)
// (3.1.5) extension in mainMotifComparison.cpp (BLiC)
// (3.2.0) implement positional prior (PSP) for mainMotifSampler.cpp 
#define VERSION "3.2.1 - version March 18 2013 by Marleen Claeys"

// type defintions for convenience
typedef double** Matrix;
typedef int** Counts;
typedef vector<int> Sequence;
typedef vector<double> ScoreVector;

// define strand modes
enum strand_modes {minus_strand=-1, both=0, plus_strand=1};

#endif

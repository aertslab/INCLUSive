#ifndef definitions_include_declared
#define definitions_include_declared

#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <cstring>

using namespace std;

#define READ 0
#define WRITE 1

#define SINGLE 0
#define BOTH 1

// type defintions for convenience
typedef double (* Matrix)[4];
typedef int  (* Counts)[4];
typedef vector<int> Sequence;
typedef vector<double> ScoreVector;

// define strand modes
enum strand_modes {minus_strand=-1, both=0, plus_strand=1};

#endif

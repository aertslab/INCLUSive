#ifndef distribution_include_declared
#define distribution_include_declared

#include "inclusive.h"
#include "RandomNumber.h"
#include <vector>

class Distribution {
 private:
  bool _normalized;
  vector<double> *_dist;
  RandomNumber _rn;
  bool _Normalize();

 public:
	// constructors
  Distribution(int length, float *pArray);
  Distribution(int length, int *pArray);
  Distribution(int length, double *pArray);
  Distribution(vector<float> *pArray);
  Distribution(vector<double> *pArray);
  Distribution(vector<int> *pArray);
  Distribution(ScoreVector *pArray, int start, int length); // 3.1.5
  // destructor
  ~Distribution();

  // inspectors
  int TakeSample();
  int SelectMax();
  bool IsNormalized(){ return _normalized; };
  vector<double> * GetVector(){return _dist;}; // 3.1.5

  // adaptors
  Distribution *ApplyMask(int start, int length);
  
};

#endif

#ifndef randomnumber_include_declared
#define randomnumber_include_declared

#include <iostream>

#define  RAND1  30269
#define  RAND2  30307
#define  RAND3  30323  
  

class RandomNumber{
 private:
  static long _lrand1;
  static long _lrand2;
  static long _lrand3;

  long _lcg(int a, int c, long M, long R);
  float _mod(float a, float b);

 public:
  RandomNumber(){};
  ~RandomNumber(){};
  static void InitRandomNumber();

  float GetUniform();
  float GetExponential(float lambda);
  
};

#endif

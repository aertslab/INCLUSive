#ifndef spacerpwm_include_defined
#define spacerpwm_include_defined

#include "inclusive.h"

class SpacerPWM{
  private:
  int _length1;
  int _length2;
  int _spacer;
  Matrix _pMatrix1;
  Matrix _pMatrix2;
  double _pPseudo[4];
  double _score;
  string * _id;
  string * _consensus;

  void _ComputeConsensus();

  public:
  // constructor
  SpacerPWM(int W1, Matrix pM1, int spacer, int W2, Matrix pM2);
  
  // destructor
  ~SpacerPWM();
  
  // inspectors
  string * GetID(){return _id;};
  string * GetConsensus(){ return _consensus;};
  int Length(){ return (_length1 + _spacer + _length2);};
  int Spacer(){ return _spacer;};
  int Length1(){ return _length1;};
  int Length2(){ return _length2;};

  // adaptors
  SpacerPWM * SetConsensus(string * cons);
  SpacerPWM * SetID(string * id);
  SpacerPWM * SetScore(double sc);
};


#endif

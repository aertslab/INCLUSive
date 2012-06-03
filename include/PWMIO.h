#ifndef pwmio_include_declared
#define pwmio_include_declared

#include "inclusive.h"
#include "PWM.h"

class PWMIO {
 private:
  bool _isOpenReading;
  bool _isOpenWriting;
  ifstream _ifs;
  ofstream _ofs;
  string _pLine;
  string * _pError; 

 public: 
    // constructor
  PWMIO(string *fileName, int type);

  // destructor
  ~PWMIO();

  // inspectors
  bool IsOpen();
  string * GetError(){return _pError;}

  PWM * ReadMatrix();
  int WriteMatrix(PWM *matrix);
  void Close();

};

#endif

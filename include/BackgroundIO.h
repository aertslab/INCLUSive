#ifndef backgroundio_include_declared
#define backgroundio_include_declared

#include "inclusive.h"
#include "BackgroundModel.h"

class BackgroundIO{
 private:
  bool _isOpenReading;
  bool _isOpenWriting;
  ifstream _ifs;
  ofstream _ofs;
  string _pLine;
  string * _pError;
  
 public:
  // constructor 
  BackgroundIO(string fileName, int type);
  // destructor
  ~BackgroundIO();

  // inspectors
  bool IsOpen();
  string * GetError(){return _pError;}
  
  BackgroundModel * ReadBackgroundModel();
  int WriteBackgroundModel(BackgroundModel *bgM);
  void Close();

};

#endif

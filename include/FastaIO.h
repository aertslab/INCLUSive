#ifndef fastaio_include_declared
#define fastaio_include_declared

#include "inclusive.h"
#include "SequenceObject.h"


class FastaIO {

 private:
  ifstream _ifs;
  bool _isOpen;
  bool _hasNext;
  string _pLine;

 public:
  FastaIO(string *fileName);
  ~FastaIO();
  
  // inspectors
  bool IsOpen();
  bool HasNext();
  SequenceObject * NextSequence();

  // adaptors
  void Close();

};

#endif

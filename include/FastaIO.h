// 21 july 2009 - 3.1.5

#ifndef fastaio_include_declared
#define fastaio_include_declared

#include "inclusive.h"
#include "SequenceObject.h"
#include <sstream> //

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
  double * ReadDirichlet(); // for MotifComparison
  string * ReadSeqID(); // for MotifSampler::UpdatePspScores
  ScoreVector * LoadPspData(int L, int w, bool skip);// for MotifSampler::UpdatePspScores

  // adaptors
  void Close();

};

#endif

#ifndef backgroundmodel_include_declared
#define backgroundmodel_include_declared

#include "inclusive.h"
#include <string>

class BackgroundModel {
 private:
  int _order;
  int _length;
  string _organism;
  string *_pFile;
  double *_snf;
  double *_oligoFrequencyMatrix;
  double (* _transitionMatrix)[4];

 public:
  // constructor
  BackgroundModel(int order, double (*pTrans)[4], double *pFreq, double *pSnf);
  // destructor
  ~BackgroundModel();
  
  // inspectors
  int GetOrder();
  string * GetOrganism();
  string * GetFileName(){return _pFile;};
  
  // double (*)[4] GetTransitionMatrix();
  double * GetSNF();
  double * GetOligoFrequencyMatrix();
  double GetSnfValueAt(int i);
  double GetTransitionMatrixValueAt(int i, int j);
  double GetOligoFrequencyValueAt(int i);

  // adaptors
  BackgroundModel * SetSequences(string * pFile){ _pFile = pFile; return this;};
  BackgroundModel * SetOrganism(string * pOrganism){  _organism = *pOrganism; return this;};
};

#endif

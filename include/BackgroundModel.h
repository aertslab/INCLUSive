#ifndef backgroundmodel_include_declared
#define backgroundmodel_include_declared

#include "inclusive.h"
#include <string>

class BackgroundModel {
 private:
  int _order;
  int _length;
  string *_organism;
  string *_pFile;
  double *_snf;
  double *_oligoFrequencyMatrix;
  double ** _transitionMatrix;

 public:
  // constructor
  BackgroundModel(int order, double **pTrans, double *pFreq, double *pSnf);
  // destructor
  ~BackgroundModel();
  
  // inspectors
  int GetOrder(){return _order;};
  string * GetOrganism(){return _organism;};
  string * GetFileName(){return _pFile;};
	double * GetSNF(){return _snf;};
  double * GetOligoFrequencyMatrix(){return _oligoFrequencyMatrix;};
  double GetSnfValueAt(int i);
  double GetTransitionMatrixValueAt(int i, int j);
  double GetMaxTransitionMatrixValueAt(int j);
  double GetMinTransitionMatrixValueAt(int j);
  double GetOligoFrequencyValueAt(int i);

  // adaptors
  BackgroundModel * SetSequences(string * pFile){ _pFile = pFile; return this;};
  BackgroundModel * SetOrganism(string * pOrganism){  _organism = pOrganism; return this;};

};

#endif

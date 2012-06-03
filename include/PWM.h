// 21 july 2009 - 3.1.5

#ifndef pwm_include_declared
#define pwm_include_declared

#include "inclusive.h"
#include "InstanceMap.h"


class PWM{

 private:
  int _length;
  Matrix _pMatrix;
  double _pPseudo[4];
  double _score;
  double _weight; // (2009/07/19 - MotifComparison BLiC)
  string * _id;
  string * _consensus;

  // (2009/09/03 - MotifComparison KL)
  int _MIshift;

  void _ComputeConsensus();

 public: 
  PWM(PWM &pCopy);
  PWM(int W, Matrix pM);
  PWM(int W, Matrix pM, double *snf);
  PWM(int W, Counts pC);
  PWM(int W, Counts pC, double *snf);
  PWM(InstanceMap * pInstances, double *snf, int w);
  ~PWM();

  Matrix GetMatrix() const { return _pMatrix;};
  int Length() const { return _length; };
  string * GetID() const {return _id;};
  double Score() const { return _score;};
  string * GetConsensus();
  void StderrPrintMatrix();
  double GetWeight(){return _weight;};//
  int GetMIshift(){return _MIshift;};// 2009/09/03

  double GetValueAt(int i, int j);
  double GetPseudoCountAt(int i);
  char GetConsensusSymbolAt(int index);
  
  double ConsensusScore();
  double InformationContent(double *snf);
  double MaxScore();

  double MutualInformation(PWM *sbjct, int shift);
  double * BLiC(PWM *dbMatrix, int shift, double *dirichlet, double *bgsnf, bool a2a1);
  double WeightedKullbackLeiberDistance(PWM * sbjct, vector<double>* pVec);

  // adaptors
  PWM * SetID(string * id);
  PWM * SetConsensus(string * cons);
  PWM * RebuildMatrix(int W, Counts pC, double *snf);
  PWM * RebuildMatrix(int W, Matrix pM, double *snf);
  PWM * RebuildMatrix(int W, Counts pC);
  PWM * RebuildMatrix(int W, Matrix pM);
  PWM * SetScore(double sc);
  void SetWeight(double w){ _weight = w; return;};//

  // sub matrix selection
  PWM * SubMatrix(int index, int length);
  PWM * ReverseSubMatrix(int index, int length);
  
  // matrix shuffling
  void ShuffleMatrix();
  PWM * NewShuffledMatrix();
  
};


#endif

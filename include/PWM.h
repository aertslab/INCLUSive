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
  string * _id;
  string * _consensus;

  void _ComputeConsensus();

 public: 
  PWM(PWM &pCopy);
  PWM(int W, Matrix pM);
  PWM(int W, Matrix pM, double *snf);
  PWM(int W, Counts pC);
  PWM(int W, Counts pC, double *snf);
  PWM(InstanceMap * pInstances, double *snf, int w);
  ~PWM();

  Matrix GetMatrix(){ return _pMatrix;};
  int Length(){ return _length; };
  string * GetID(){return _id;};
  double Score(){ return _score;};
  string * GetConsensus();
  void StderrPrintMatrix();

  double GetValueAt(int i, int j);
  double GetPseudoCountAt(int i);
  char GetConsensusSymbolAt(int index);
  
  double ConsensusScore();
  double InformationContent(double *snf);
  double MaxScore();

  double MutualInformation(PWM *sbjct, int shift);
  double WeightedKullbackLeiberDistance(PWM * sbjct, vector<double>* pVec);

  // adaptors
  PWM * SetID(string * id);
  PWM * SetConsensus(string * cons);
  PWM * RebuildMatrix(int W, Counts pC, double *snf);
  PWM * RebuildMatrix(int W, Matrix pM, double *snf);
  PWM * RebuildMatrix(int W, Counts pC);
  PWM * RebuildMatrix(int W, Matrix pM);
  PWM * SetScore(double sc);

  // sub matrix selection
  PWM * SubMatrix(int index, int length);
  
};


#endif

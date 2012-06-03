#ifndef utilities_include_defined
#define utilities_include_defined

#include "inclusive.h"
#include "SequenceObject.h"
#include "PWM.h"
#include "BackgroundModel.h"

/*
 * Here the INCLUSive namespace is defined where all 
 * auxillary function are declared that could not be
 * assigned to a specific class 
 */

namespace INCLUSIVE {

  // general mathematical functions  
  double SumArray(const double *pt, const int start, const int number);
  double SumArray(vector<double>::iterator pt, const int start, const int number);
  double SumArray(vector<double> *pt, const int start, const int number);
  int MaxIndex(vector<double> *pt);
  double LogGamma(const double x);
  double LogNormDirichlet(vector<double>* pVec);
  double LogDirichlet(vector<double>* pParam, vector<double>* pVec);

  // utilities to transform short sequences to indices and vice versa
  int Site2Index(SequenceObject *pSeq, strand_modes s, int start, int order);
  void Index2Site(string &rStr, int i, int L);
  
	// specific sequence scoring functions 
  void SegmentLogMatrixScore(ScoreVector &pmx, SequenceObject *pSeq, strand_modes s, PWM *pMatrix);
  void SegmentLogBackgroundScore(ScoreVector &pbx, SequenceObject *pSeq, strand_modes s, BackgroundModel *pBgM, int w);
  double SequenceLogBackgroundScore(SequenceObject *pSeq, strand_modes s, BackgroundModel *pBgM);
  double SequenceLogBackgroundScore(SequenceObject *pSeq, int start, int length, strand_modes s, BackgroundModel *pBgM);

  double ComputeNInstancesProbability(ScoreVector *pExpWx, int n, int L, int W);
  void ComputePriorDistribution(ScoreVector &pPrior, int n, double prior);

  // utility to reset scoring vector
  void ResetScoreVector(ScoreVector &pVec);
  void ResetScoreVector(ScoreVector &pVec, int value);

  // faculteit
  int fac(int n);
  
}

#endif

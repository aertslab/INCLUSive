#ifndef sequencecomputation_include_declared
#define sequencecomputation_include_declared

#include "SequenceObject.h"
#include "MaskVector.h"
#include "inclusive.h"
#include "utilities.h"
#include "Distribution.h"
#include <cmath>

class SequenceComputation{
  private:
    SequenceObject *_pSeqObj; 
    int _wlength, _length;
    bool correct;
  
    // internal score vectors
    ScoreVector *_pMotifScore;
    double _dLogP0;
    ScoreVector *_pBackgroundScore;
    ScoreVector *_pExpWx;
    MaskVector *_pMask;
  
    // reverse strand score vectors
    ScoreVector *_pRevMotifScore;
    double _dRevLogP0;
    ScoreVector *_pRevBackgroundScore;
    ScoreVector *_pRevExpWx;
    MaskVector *_pRevMask;
 
    ScoreVector *_pCopyProbDistr;
    ScoreVector *_pRevCopyProbDistr;
    ScoreVector *_pPriorDistr;
  
  public:
    // constructor
    SequenceComputation(SequenceObject *pSeq);
    
    // destructor
    ~SequenceComputation();
  
    // update motif length
    void SetMotifLength(int wLength);
  
    // update scores
    void UpdateInstanceMotifScore(PWM *pMotif, strand_modes STRAND);
    void UpdateInstanceBackgroundScore(BackgroundModel *pBg, int wLength, strand_modes STRAND);
    void UpdateSequenceBackgroundScore(BackgroundModel *pBg, strand_modes STRAND);
    void UpdateInstanceExpWx(strand_modes STRAND); 
  
    // mask operations
    void UpdateMask(int start, int w, int value, strand_modes STRAND);
    void SetMask(MaskVector *pMask, strand_modes STRAND);
    void ApplyMask(strand_modes STRAND);
    
    void UpdateCopyProbability(double prior, strand_modes STRAND);
    void UpdateFixedSizeCopyProbability(int maxN, double prior, strand_modes STRAND);
    void FixCopyProbability(ScoreVector *pCopyProbValues, strand_modes STRAND);
  
    // start positions of selected instances
    void SampleInstanceStart(vector<int> & pAlignmentVector, int n, strand_modes STRAND);
    void SampleUniformInstanceStart(vector<int> & pAlignmentVector, int n, strand_modes STRAND);
    void SelectBestInstanceStart(vector<int> & pAlignmentVector, int n, strand_modes STRAND);

    // inspectors
    int GetEstimatedNumberInstances(strand_modes STRAND);
    double GetWxAt(int ndx, strand_modes STRAND);
    MaskVector * GetMask(strand_modes STRAND);    
    ScoreVector * GetCopyProbability(strand_modes STRAND);
    double GetCopyProbabilityAt(int nbr, strand_modes STRAND);
    double GetLogBackgroundScore(strand_modes STRAND);
    SequenceObject * ParentSequence(){return _pSeqObj;};
    double LogLikelihoodScore(vector<int> *pAlign, strand_modes STRAND);
  
};

#endif

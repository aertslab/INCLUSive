#ifndef priormotifsamplerrun_include_declared
#define priormotifsamplerrun_include_declared

#include "MotifSamplerRun.h"
#include "inclusive.h"
#include "utilities.h"


class PriorMotifSamplerRun : public MotifSamplerRun 
{
  private:
    Matrix _pPriorMotifModel;
    vector<double>* _pMixtureCoefficients;
  
  public:
    // constructor calls base class
    PriorMotifSamplerRun(string* pFasta, strand_modes strand, PWM* pPriorMotif);

    // destructor
    ~PriorMotifSamplerRun();

    // reset the motif models
    void SetNewPriorMotifModel(PWM* pPriorMotif);

    // methods to compute the pseudo counts
    void UpdateWeightingFactor(double w){_scale = w; return;};
    void UpdateMixtureCoefficients();
    void ResetMixtureCoefficients();
    
    // methods to update motif models
    void BuildMotifFromInstanceMap();
    void BuildMotifFromReducedInstanceMap(SequenceObject* pSeq);
    
    // running procedures
    void InitializationStep(int iterations);
    void CoreSamplingStep(int iterations);
    void ConvergenceStep(int iterations);
    
    // printing
    void StderrPrintMixtureCoefficients();
    
    
    vector<double>* MixtureCoefficients(){return _pMixtureCoefficients;};
    
};


#endif

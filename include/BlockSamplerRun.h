#ifndef blocksamplerrun_include_declared
#define blocksamplerrun_include_declared

#include "MotifSamplerRun.h"
#include "inclusive.h"
#include "utilities.h"

class BlockSamplerRun : public MotifSamplerRun 
{
  private:
    map<string, double> *_pWeightMap;
    SequenceObject* _pRootSequence;
    
  public:
    BlockSamplerRun(string* pFasta, string* pRootID, strand_modes strand);
    ~BlockSamplerRun();
  
    void SetSequenceWeight();
    
    void InitializationStep(int iterations);
    void CoreSamplingStep(int iterations);
    void ConvergenceStep(int iterations);
  
};

#endif

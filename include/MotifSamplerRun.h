#ifndef motifsamplerrun_include_declared
#define motifsamplerrun_include_declared

// #include <fstream>
#include "utilities.h"
#include "inclusive.h"
#include "SequenceObject.h"
#include "SequenceComputation.h"
#include "PWM.h"
#include "BackgroundModel.h"
#include "Instance.h"
#include "InstanceMap.h"
#include "Distribution.h"
#include "GFFWriter.h"

#include <list>

// some useful type definitions 
typedef list<SequenceObject*> SequenceSet;
typedef list<SequenceObject *>::iterator SeqIterator;
typedef map<SequenceObject *, SequenceComputation *>::iterator MapIterator;
typedef map<SequenceObject *, SequenceComputation *>::value_type MapValue;


class MotifSamplerRun{

 protected:
  // sequence data
  // store sequences in list
  SequenceSet *_pSequenceList;
  
  // store computation objects in map, where the key is a pointer to the sequence object 
  // and the associated value points to the computation object 
  map<SequenceObject *, SequenceComputation *> * _pComputationMap;
  strand_modes _strand;
  int _nbrSequences;

  // instances
  InstanceMap * _pInstanceMap;
  int _nbrInstances;
  int _nbrSequencesWithInstances;

  // parameters
  int _w;
  double _prior;
  int _overlap;
  string *_pId;
  double _scale;

 // motif model
  double* _pPseudoCounts;
  Counts _pLocalCounts;
  PWM *_localMotif;

  // background model
  BackgroundModel *_pBgModel;
  
  // local methods
  void _LoadSequences(string * pFastaFile, strand_modes strand);
  void _ClearInstanceMap();

 public:

  // constructor
  MotifSamplerRun(string * pFastaFile, strand_modes strand);
  
  // destructor
  ~MotifSamplerRun();
  
  // methods to set parameters
  void SetMotifLength(int wLength); // setting the motif length will change the masks
  void SetInstancePrior(int prior){_prior = prior; return;};
  void SetOneInstancePrior(double prior){_prior = prior; return;};
  void SetInstanceID(string *pSiteId){_pId = pSiteId; return;};
  void SetOverlap(int overlap){_overlap = overlap; return;};
  void SetBackgroundModel(BackgroundModel *pBg);
  void UpdatePseudoCounts(double* psCounts);
  void UpdatePseudoCounts(double* psCounts, double scale);
  
  // instance map creation/manipulation
  void InitFixedSizeInstanceMap(int nbr);
  void SampleFixedSizeInstanceMap(int nbr);
  void SampleMaxSizeInstanceMap(int nbr);
  void SampleInstanceMap();
  void SelectFixedSizeBestInstanceMap(int nbr);
  void SelectMaxSizeBestInstanceMap(int nbr);
  void SelectBestInstanceMap();
  void ShiftInstanceMap(int maxShift);

  // motif scores
  void UpdateMotifScores();
  void UpdateMotifScoresFromReducedInstanceMap();
  void UpdateCopyProbabilityDistribution();
  
  // background model scores
  void UpdateBackgroundScores();
  bool UpdateBackgroundScore(string* pID, BackgroundModel *pBgM);

  // motif manipulation -> stored in _localMotif
  void BuildMotifFromInstanceMap();
  void BuildMotifFromReducedInstanceMap(SequenceObject *pSeq);
  PWM * GetMotifModel();
  double LogLikelihoodScore();

  // list inspectors
  int NumberOfSequences(){ return _nbrSequences;};
  int NumberOfInstances(){ return _nbrInstances;};
  int NumberOfSequencesWithInstances(){ return _nbrSequencesWithInstances;};

  // different steps in procedure
  void InitializationStep(int iterations);
  void CoreSamplingStep(int iterations, int shiftTime, int maxShift);
  void CoreMaxSizeSamplingStep(int maxNbr, int iterations, int shiftTime, int maxShift);
  void ConvergenceStep(int iterations);
  void MaxSizeConvergenceStep(int maxNbr, int iterations);
  
  // masking vectors
  void ResetMasks();
  void UpdateMasksFromInstanceMap();

  // reporting and visualization
  void StderrPrintMotifInfo();
  void PrintInstanceMap(GFFWriter *pGFFout, string* pSource);
  void StderrPrintInstanceMap();

};


#endif
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
  //double _prior;
  vector<Distribution*> * _pPriorDistributions; 
  int _overlap;
  string *_pId;
  double _scale;

 // motif model
  double* _pPseudoCounts;
  double** _pLocalCounts; // double instead of Counts because of PSP
  PWM *_localMotif;

  // background model
  BackgroundModel *_pBgModel;
  
  // local methods
  void _LoadSequences(string * pFastaFile, strand_modes strand);
  void _ClearInstanceMap();

  // output
  GFFWriter * _pgffio;

 public:

  // constructor
  MotifSamplerRun(string * pFastaFile, strand_modes strand, GFFWriter * pgff);
  
  // destructor
  ~MotifSamplerRun();
  
  // methods to set parameters
  void SetMotifLength(int wLength); // setting the motif length will change the masks
  void SetInstanceID(string *pSiteId){_pId = pSiteId; return;};
  void SetOverlap(int overlap){_overlap = overlap; return;};
  void SetBackgroundModel(BackgroundModel *pBg);
  void UpdatePseudoCounts(double* psCounts);
  void UpdatePseudoCounts(double* psCounts, double scale);
  bool LoadNbrInstInfo(int maxM, string * priorinput, bool sample); // 
  bool LoadPspScores(string * pPspFile, bool bPSPs); //
  
  // instance map creation/manipulation
  void InitFixedSizeInstanceMap(int nbr);
  void SampleFixedSizeInstanceMap(int nbr);
  void SampleInstanceMap();
  void SelectFixedSizeBestInstanceMap(int nbr);
  void SelectBestInstanceMap();
  void ShiftInstanceMap(int maxShift);
  void ExtendLeftInstanceMap(int pos);
  void ExtendRightInstanceMap(int pos);
  int CheckLeftExtension(double threshold);
  int CheckRightExtension(double threshold);
  int CheckLeftExtension(double threshold, int window);
  int CheckRightExtension(double threshold, int window);
  
  // motif scores
  void UpdateMotifScores();
  void UpdateMotifScoresFromReducedInstanceMap();
  void UpdateCopyProbabilityDistribution();
  
  // background model scores
  void UpdateBackgroundScores();
  bool UpdateBackgroundScore(const string& pID, BackgroundModel *pBgM);

  // motif manipulation -> stored in _localMotif
  void BuildMotifFromInstanceMap();
  void BuildMotifFromReducedInstanceMap(SequenceObject *pSeq);
  PWM * GetMotifModel();
  double LogLikelihoodScore();

  // list inspectors
  int NumberOfSequences(){ return _nbrSequences;};
  int NumberOfInstances(){ return _nbrInstances;};
  int NumberOfSequencesWithInstances(){ return _nbrSequencesWithInstances;};
  SequenceComputation * FindSequence(string * id);

  // different steps in procedure
  void InitializationStep(int iterations);
  void CoreSamplingStep(int iterations, int shiftTime, int maxShift);
  void ConvergenceStep(int iterations);
  
  // masking vectors
  void ResetMasks();
  void ResetMasks(int wNew);
  void UpdateMasksFromInstanceMap();

  // reporting and visualization
  void StderrPrintMotifInfo();
  void PrintInstanceMap(GFFWriter *pGFFout, string* pSource);
  void StderrPrintInstanceMap();
};


#endif

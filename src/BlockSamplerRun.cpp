#include "BlockSamplerRun.h" 
#include "inclusive.h"
#include "utilities.h"

using namespace INCLUSIVE;


BlockSamplerRun::BlockSamplerRun(string* pFasta, string* pRootID, strand_modes strand)
  : MotifSamplerRun(pFasta, strand, NULL)
{
  
  // set the root sequence pointer
  MapIterator mi = _pComputationMap->begin();
  SequenceObject* pSeq = NULL;
  while ( mi != _pComputationMap->end() )
  {
    pSeq = (*mi).first;
    if ( *(pSeq->GetID()) == *pRootID )
    {
      _pRootSequence = pSeq;
      _pRootComputation = (*mi).second; // 
      break;
    }
    mi++;
  }
  pSeq = NULL;
  cerr << "DEBUG:  Root sequence: " << *(_pRootSequence->GetID()) << endl;
}

BlockSamplerRun::~BlockSamplerRun()
{
  _pRootSequence =  NULL; 
  _pRootComputation = NULL;
}

// extra function (3.1.5), this only needs to be done once, not every iteration
void
BlockSamplerRun::FixNbrInstForRoot()
{
    // hard code the copy probability to always select 1 instance
  ScoreVector* pFixCopy = new ScoreVector(2,0);
  (*pFixCopy)[0] = 0.001;
  (*pFixCopy)[1] = 0.999;
  // a way to use the same function as MotifSamplerRun
  _pRootComputation->LoadNbrInstFixed(pFixCopy);
  
  // cleanup
  delete pFixCopy; pFixCopy = NULL;
 
  return; 
}

void 
BlockSamplerRun::InitializationStep(int iterations)
{
  // check the number of iterations
  MapIterator mi = _pComputationMap->begin();
  SequenceComputation *pComp = NULL;

  cerr << "Initialize Instance Map" << endl;
  InitFixedSizeInstanceMap(1);
  StderrPrintInstanceMap();
  
  cerr << "Starting initial 1 copy step" << endl;
  for (int i = 0; i < iterations; i++)
  {
    // loop over all sequences and update scores
    mi = _pComputationMap->begin();
    while (mi != _pComputationMap->end())
    {
      pComp = (*mi).second;

      // build motif model from instances in instance map excluding current sequence 
      BuildMotifFromReducedInstanceMap((*mi).first);

      // update motif scores
      pComp->UpdateInstanceMotifScore(_localMotif, _strand);
      pComp->UpdateInstanceExpWx(_strand);

      // next sequence
      mi++;
    }

    // sample exactly one instance in each sequence 
    SampleFixedSizeInstanceMap(1);

    // print information on screen
    cerr << i << "\t";
    StderrPrintMotifInfo();

  }

  pComp = NULL;  
  return;
}


void
  BlockSamplerRun::CoreSamplingStep(int iterations)
{
  MapIterator mi = _pComputationMap->begin();
  SequenceComputation* pComp = NULL;
  //SequenceObject* pSeq = NULL;
  
/*
  // hard code the copy probability to always select 1 instance
  ScoreVector* pFixCopy = new ScoreVector(2,0);
  (*pFixCopy)[0] = 0.001;
  (*pFixCopy)[1] = 0.999;  
*/
  for (int i = 0; i < iterations; i++)
  {

    // loop over all sequences and update scores
    mi = _pComputationMap->begin();
    while (mi != _pComputationMap->end())
    {
      pComp = (*mi).second;
      //pSeq = (*mi).first;
      
      // build motif model from instances in instance map excluding current sequence 
      BuildMotifFromReducedInstanceMap((*mi).first);

      // update motif scores
      pComp->UpdateInstanceMotifScore(_localMotif, _strand);
      pComp->UpdateInstanceExpWx(_strand);

      // compute distribution to estimate number of instances
      pComp->UpdateCopyProbability(_strand); // 3.1.5, covers what follows
  
/*     // if this is the root sequence set the hard coded copy probability
      if ( pSeq == _pRootSequence )
      {
        pComp->FixCopyProbability(pFixCopy, _strand);
      }
      else
      {
        pComp->UpdateFixedSizeCopyProbability(1,_prior, _strand);
      }
 */
      // next sequence
      mi++;
    }

    // sample a new instance map
    //SampleMaxSizeInstanceMap(1);
    SampleInstanceMap();//will automatically be 1 (as -M=1 main)

    // check number of instances
    if (_nbrInstances < 2 || _nbrSequencesWithInstances < 2)
    {
      cerr << "--Warning-- not enough instances to proceed procedure." << endl;
      pComp = NULL;
      break;
    }

    // print information on current motif
    cerr << i << "\t";
    StderrPrintMotifInfo();

    pComp = NULL;
  }

  // delete local variables
  //delete pFixCopy;
  
  return;
}


/****************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
void
  BlockSamplerRun::ConvergenceStep(int iterations)
{
  // cerr << "DEBUG: MotifSamplerRun::ConvergenceStep" << endl;
  MapIterator mi = _pComputationMap->begin();
  //SequenceObject* pSeq = NULL;
  SequenceComputation* pComp = NULL;
/*  // hard code the copy probability to always select 1 instance
  ScoreVector * pFixCopy = new ScoreVector(2);
  (*pFixCopy)[0] = 0.001;
  (*pFixCopy)[1] = 0.999;  
*/
  for (int i = 0; i < iterations; i++)
  {
    // build motif model from all instances
    BuildMotifFromInstanceMap();

    // loop over all sequences
    mi = _pComputationMap->begin();
    while (mi != _pComputationMap->end())
    {
      //pSeq = (*mi).first;
      pComp = (*mi).second;

      // update motif scores
      pComp->UpdateInstanceMotifScore(_localMotif, _strand);
      pComp->UpdateInstanceExpWx(_strand);

      // compute distribution to estimate number of instances
      pComp->UpdateCopyProbability(_strand); // 3.1.5, covers what follows
/*
      // if this is the root sequence set the hard coded copy probability
      if (  pSeq == _pRootSequence )
      {
        pComp->FixCopyProbability(pFixCopy, _strand);
      }
      else
      {
        pComp->UpdateFixedSizeCopyProbability(1,_prior, _strand);
      }
*/
      // next sequence
      mi++;
    }

    // create new instance map from curent scores
    SelectBestInstanceMap();

    // check number of instances
    if (_nbrInstances < 2 || _nbrSequencesWithInstances < 2)
    {
      cerr << "--Warning-- not enough instances to proceed procedure." <<
        endl;
      //pSeq = NULL;
      pComp = NULL;
      break;
    }

    cerr << i << "\t";
    StderrPrintMotifInfo();

    //pSeq = NULL;
    pComp = NULL;
  }

  // final motif model
  BuildMotifFromInstanceMap();

  // clear pFixCopy
  //delete pFixCopy;
  
  return;
}

#include "SequenceComputation.h"
#include "utilities.h"
#include <cmath>

using namespace INCLUSIVE;

/******************************************************************************
  Method:       new
  Class:        SequenceComputation
  Arguments:    SequenceObject * pSeqObj
  
  Description:  creates a new computation object
  
  Date:         2012/09/19
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
                Revised : Marleen Caeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
SequenceComputation::SequenceComputation(SequenceObject * pSeqObj)
{
  // start with setting the link to parental sequence object
  _pSeqObj = pSeqObj;

  // sequence length
  _length = _pSeqObj->Length();

  // initialize vectors
  _pMotifScore = new ScoreVector(_length, 0);
  _pBackgroundScore = new ScoreVector(_length, 0);
  _pExpWx = new ScoreVector(_length, 0);
  _pMask = new MaskVector(_length);
  _pPspScore = new ScoreVector(_length, 1);

  _pRevMotifScore = new ScoreVector(_length);
  _pRevBackgroundScore = new ScoreVector(_length, 0);
  _pRevExpWx = new ScoreVector(_length, 0);
  _pRevMask = new MaskVector(_length);
  _pRevPspScore = new ScoreVector(_length, 1);

  // set initial length of copy probability to 10
  _pCopyProbDistr = new ScoreVector(10, 0);
  _pRevCopyProbDistr = new ScoreVector(10, 0);
  //_pPriorDistr = new ScoreVector(10, 0);
  _priorDistrs = NULL; // will be created in LinkNbrInstPrior
  _plogNormFactors = NULL; // will be updated based on PSP

  // initial motif length is zero
  _wlength = 0;

  // initially the format is correct
  _correct = true;

}

/******************************************************************************
  Method:       delete
  Class:        SequenceComputation
  Arguments:    none
  
  Description:  destructor
  
  Date:         2012/09/19
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
                Revised : Marleen Caeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
SequenceComputation::~SequenceComputation()
{
  // cleanup internal vectors

  // score of the plus strand
  if (_pMotifScore != NULL)
    delete _pMotifScore;
  _pMotifScore = NULL;

  if (_pBackgroundScore != NULL)
    delete _pBackgroundScore;
  _pBackgroundScore = NULL;

  if (_pExpWx != NULL)
    delete _pExpWx;
  _pExpWx = NULL;

  if (_pCopyProbDistr != NULL)
    delete _pCopyProbDistr;
  _pCopyProbDistr = NULL;

  if (_pMask != NULL)
    delete _pMask;
  _pMask = 0;
	
  if (_pPspScore != NULL)
    delete _pPspScore;
  _pPspScore = NULL;

  // scores of the minus strand
  if (_pRevMotifScore != NULL)
    delete _pRevMotifScore;
  _pRevMotifScore = NULL;

  if (_pRevBackgroundScore != NULL)
    delete _pRevBackgroundScore;
  _pRevBackgroundScore = NULL;

  if (_pRevCopyProbDistr != NULL)
    delete _pRevCopyProbDistr;
  _pRevCopyProbDistr = NULL;

  if (_pRevExpWx != NULL)
    delete _pRevExpWx;
  _pRevExpWx = NULL;

  if (_pRevMask != NULL)
    delete _pRevMask;
  _pRevMask = 0;

  if (_pRevPspScore != NULL)
    delete _pRevPspScore;
  _pRevPspScore = NULL;

  /*/ other 
  if (_pPriorDistr != NULL)
    delete _pPriorDistr;
  _pPriorDistr = NULL;*/
  //_priorDistrs = NULL; // do not delete, is done in main
  // delete, new copies were created
  if (_priorDistrs != NULL)
  { // cleanup distributions
    for(int i = 0; i < (int)_priorDistrs->size(); i++)
    { if ( (*_priorDistrs)[i] != NULL)
        delete (*_priorDistrs)[i];
    }
    delete _priorDistrs;
  }
  _priorDistrs = NULL;

  if (_plogNormFactors != NULL)
    delete _plogNormFactors;
  _plogNormFactors = NULL;

}


/******************************************************************************
  Method:       UpdateInstanceMotifScore 
  Class:        SequenceComputation
  Arguments:    PWM * pMotif, strand_modes STRAND
  
  Description:  Update the motif score of all instances on STRAND with 
								motif model pMotif
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::UpdateInstanceMotifScore(PWM * pMotif,
                                              strand_modes STRAND)
{
  // update internal motif length
  _wlength = pMotif->Length();

  // score plus strand
  if (STRAND == both || STRAND == plus_strand)
  {
    SegmentLogMatrixScore(*_pMotifScore, _pSeqObj, plus_strand, pMotif);
  }

  // score minus strand
  if (STRAND == both || STRAND == minus_strand)
  {
    SegmentLogMatrixScore(*_pRevMotifScore, _pSeqObj, minus_strand, pMotif);
  }

  return;
}


/******************************************************************************
  Method:       UpdateSequenceBackgroundScore
  Class:        SequenceComputation
  Arguments:    BackgroundModel * pBg, strand_modes STRAND
  
  Description:  Compute the background score of the full sequence
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::UpdateSequenceBackgroundScore(BackgroundModel * pBg,
                                                   strand_modes STRAND)
{
  // score plus strand
  if (STRAND == both || STRAND == plus_strand)
  {
    _dLogP0 = SequenceLogBackgroundScore(_pSeqObj, plus_strand, pBg);
  }

  // score minus strand
  if (STRAND == both || STRAND == minus_strand)
  {
    _dRevLogP0 = SequenceLogBackgroundScore(_pSeqObj, minus_strand, pBg);
  }

  return;
}


/******************************************************************************
  Method:       UpdateInstanceBackgroundScore
  Class:        SequenceComputation
  Arguments:    BackgroundModel * pBg, int wLength, strand_modes STRAND
  
  Description:  Score all instances of length wLength on STRAND with the 
                background model pBg

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::UpdateInstanceBackgroundScore(BackgroundModel * pBg,
                                                   int wLength,
                                                   strand_modes STRAND)
{
  // score plus strand
  if (STRAND == both || STRAND == plus_strand)
  {
    SegmentLogBackgroundScore(*_pBackgroundScore, _pSeqObj, plus_strand, pBg,
                              wLength);
  }

  // score minus strand
  if (STRAND == both || STRAND == minus_strand)
  {
    SegmentLogBackgroundScore(*_pRevBackgroundScore, _pSeqObj, minus_strand,
                              pBg, wLength);
  }

  return;
}


/******************************************************************************
  Method:       SetMask
  Class:        SequenceComputation
  Arguments:    MaskVector * pMask, strand_modes STRAND
  
  Description:  define the internal _pMask with an external pointer
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::SetMask(MaskVector * pMask, strand_modes STRAND)
{
  if (pMask != NULL && STRAND == plus_strand)
  {
    _pMask = pMask;
  }
  else if (pMask != NULL && STRAND == minus_strand)
  {
    _pRevMask = pMask;
  }
  return;
}

/******************************************************************************
  Description:  LinkNbrInstPrior : link pointer and set _maxPriorSize
	  Author_Date:  MC_2013/03/11 : create local _priorDistr
******************************************************************************/
void 
SequenceComputation::LinkNbrInstPrior(int max, vector<Distribution*>* distrs)
{ 
  // set the maximum searchspace
  _maxInst = max;
  // in case f : define _p(Rev)CopyProbDistr (only once)
  // copy needed as format from Distribution to ScoreVector
  if ((*distrs)[0] != NULL)
  {
    ScoreVector *fixedDist = (*distrs)[0]->GetVector();
    for(int i = 0; i < (int)fixedDist->size(); i++)
    {
      if (i <= 9) // replace (see constructor new ScoreVector(10, 0);)
        (*_pCopyProbDistr)[i] = (*fixedDist)[i];
      else // add
        _pCopyProbDistr->push_back((*fixedDist)[i]);
    }
    for(int i = 0; i < (int)fixedDist->size(); i++)
    {
      if (i <= 9) // replace (see constructor new ScoreVector(10, 0);)
        (*_pRevCopyProbDistr)[i] = (*fixedDist)[i];
      else // add
        _pRevCopyProbDistr->push_back((*fixedDist)[i]);
    }
    _priorDistrs = NULL; // do not link
  }
  else // other cases where we need the prior information
  {
    // link the pointer 
    //_priorDistrs = distrs;
    // create/store a copy of the global variable into a local vector
    // (because its contents may be changed sequence-specific by PSP)
    _priorDistrs = new vector<Distribution *>;
    for(int i = 0; i < (int)distrs->size(); i++)
    { // create a copy of the Distribution
      Distribution * pPriorDist = NULL;
      if ((*distrs)[i] != NULL) 
        pPriorDist = new Distribution((*distrs)[i]->GetVector());
      _priorDistrs->push_back(pPriorDist);
      pPriorDist = NULL;
    }
    // create the vector to store NORM-factors (computed based on PSP)
    _plogNormFactors = new ScoreVector(_maxInst+1,0); 
    (*_plogNormFactors)[0] = 0; // reset, entry is never used
    // set the default NORM-factors to "C" (number of possible combinations)
      // note : C code comes from utilities->ComputeNInstancesProbability)
    // this C is used when PSP-is not defined (comes back to original design)
    double logC;
    for(int n = 1; n < _maxInst+1; n++)
    {
      logC = log((double) (_length - (n * _wlength) +1));
      for (int i = 2; i <= n; i++)
        logC += log((double) (_length - (n * _wlength) + i)) - log((double) i);
      (*_plogNormFactors)[n] = logC;
    }
/*
	  cerr << "debug-followup-SequenceComputation::LinkNbrInstPrior:" << 
		  "_plogNormFactors default =(";
	  for (int i = 0; i < _maxInst+1; i++)
		  cerr << (*_plogNormFactors)[i] << " ";
	  cerr << ")" << endl;
	  */
  }
  //cerr << "debug--SequenceComputation::LinkNbrInstPrior - end" << endl;
  return;
}


/******************************************************************************
  Method:       UpdateInstanceExpWx
  Class:        SequenceComputation
  Arguments:    strand_modes STRAND
  
  Description:  set the score of all instances based on the resp. motif scores
                and the background scores
  
  Date:         2012/09/19
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
	            Revised: Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
void
SequenceComputation::UpdateInstanceExpWx(strand_modes STRAND, bool psp)
{
  // define some iterators
  vector < double >::iterator iter1;
  vector < double >::iterator iter2;
  vector < double >::iterator iter3;
  vector < int >::iterator iter4;
  vector < double >::iterator iter5; // to iterate through PSP

  // plus strand
  if (STRAND == both || STRAND == plus_strand)
  {
    iter1 = _pMotifScore->begin();
    iter2 = _pBackgroundScore->begin();
    iter3 = _pExpWx->begin();
    iter4 = _pMask->begin();
    iter5 = _pPspScore->begin();

    while (iter1 != _pMotifScore->end())
    {
      if (psp)
      { (*iter3) = exp((*iter1) - (*iter2)) * (*iter4)* (*iter5);}
      else
      { (*iter3) = exp((*iter1) - (*iter2)) * (*iter4);}
      // cerr << (*iter3) << " ";
      // augment iterators
      iter1++;
      iter2++;
      iter3++;
      iter4++;
      iter5++;
    }
    // cerr << endl;
  }
  // minus strand
  if (STRAND == both || STRAND == minus_strand)
  {
    iter1 = _pRevMotifScore->begin();
    iter2 = _pRevBackgroundScore->begin();
    iter3 = _pRevExpWx->begin();
    iter4 = _pRevMask->begin();
    iter5 = _pRevPspScore->begin();

    while (iter1 != _pRevMotifScore->end())
    {
      if (psp)
      { (*iter3) = exp((*iter1) - (*iter2)) * (*iter4)* (*iter5);}
      else
      { (*iter3) = exp((*iter1) - (*iter2)) * (*iter4);}
      // augment iterators
      iter1++;
      iter2++;
      iter3++;
      iter4++;
      iter5++;
    }

  }
  return;
}


/******************************************************************************
  Method:       UpdateCopyProbability
  Class:        SequenceComputation
  Arguments:    double prior, strand_modes STRAND
  
  Description:  compute the complete probability distribution to estimate
                the number of motif instances
								
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
/*void
SequenceComputation::UpdateCopyProbability(double prior, strand_modes STRAND)
{
  // reset copy probability arrays
  vector < double >::iterator iter;
  for (iter = _pCopyProbDistr->begin(); iter != _pCopyProbDistr->end(); ++iter)
  {
    (*iter) = 0;
  }
  for (iter = _pRevCopyProbDistr->begin(); iter != _pRevCopyProbDistr->end(); ++iter)
  {
    (*iter) = 0;
  }
	
  // reset prior distribution
  for (iter = _pPriorDistr->begin(); iter != _pPriorDistr->end(); ++iter)
  {
    (*iter) = 0;
  }
	
  // initialize with 1-prior and prior
  (*_pPriorDistr)[0] = 1 - prior;
  (*_pPriorDistr)[1] = prior;

	// some variables
  double nextValue, lastValue = 1;
  uint copyNbr = 1;
  vector < double >*pPs = new vector < double >;
  pPs->push_back(_dLogP0);

  // positive strand
  if (STRAND == both || STRAND == plus_strand)
  {
    while (_length > (int) copyNbr * _wlength && lastValue >= 0.001)
    {
      // compute next value of CopyProbability
      nextValue =
        _dLogP0 + ComputeNInstancesProbability(_pExpWx, copyNbr, _length,
                                               _wlength);

      //if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      if (finite(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      {
        // add element to pPs
        pPs->push_back(nextValue);

        // intermediate values
        double *pps = new double[copyNbr + 1];
        double stot = 0;
        for (uint j = 0; j < (copyNbr + 1); j++)
        {
          pps[j] = (*pPs)[j] + log((*_pPriorDistr)[j]);
          stot += exp(pps[j] - _dLogP0);
        }
        stot = _dLogP0 + log(stot);

        // update copy probability distribution
        for (uint j = 0; j < copyNbr; j++)
        {
          (*_pCopyProbDistr)[j] = exp(pps[j] - stot);
        }
        // add new entry to copy probability distribution
        if (copyNbr < _pCopyProbDistr->size() - 1)
        {
          (*_pCopyProbDistr)[copyNbr] = exp(pps[copyNbr] - stot);
        }
        else
        {
          _pCopyProbDistr->push_back(exp(pps[copyNbr] - stot));
        }

        // define last value in probability distribution
        lastValue = (*_pCopyProbDistr)[copyNbr];

        // delete intermediate variables
        delete[] pps;
        pps = NULL;

        // update prior by adding next value to prior array
        if (copyNbr < _pPriorDistr->size() - 1)
        {
          (*_pPriorDistr)[copyNbr + 1] = (*_pPriorDistr)[copyNbr] / 4;
        }
        else
        {
          _pPriorDistr->push_back(((*_pPriorDistr)[copyNbr] / 4));
        }
				
				// normalisation
        stot = 0;
        for (uint j = 0; j < _pPriorDistr->size(); j++)
          stot += (*_pPriorDistr)[j];

        for (uint j = 0; j < _pPriorDistr->size(); j++)
          (*_pPriorDistr)[j] = (*_pPriorDistr)[j] / stot;

      }
      else
      {
        _pCopyProbDistr->push_back(0);
        lastValue = 0;
      }

      // move to next
      copyNbr++;
    }
  }

  // negative strand
  if (STRAND == both || STRAND == minus_strand)
  {
    copyNbr = 1;
    lastValue = 1;
    vector < double >*pRevPs = new vector < double >;
    pRevPs->push_back(_dLogP0);

    while (_length > (int) copyNbr * _wlength && lastValue >= 0.001)
    {
      // compute next value of CopyProbability
      nextValue =
        _dLogP0 + ComputeNInstancesProbability(_pRevExpWx, copyNbr, _length,
                                               _wlength);

      // if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      if (finite(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      {
        // add element to pRevPs
        pRevPs->push_back(nextValue);

        // intermediate values
        double *pps = new double[copyNbr + 1];
        double stot = 0;
        for (uint j = 0; j < (copyNbr + 1); j++)
        {
          pps[j] = (*pRevPs)[j] + log((*_pPriorDistr)[j]);
          stot += exp(pps[j] - _dLogP0);
        }
        stot = _dLogP0 + log(stot);

        // update copy probability distribution
        for (uint j = 0; j < copyNbr; j++)
        {
          (*_pRevCopyProbDistr)[j] = exp(pps[j] - stot);
        }
        // add new entry to copy probability distribution
        if (copyNbr < _pRevCopyProbDistr->size() - 1)
        {
          (*_pRevCopyProbDistr)[copyNbr] = exp(pps[copyNbr] - stot);
        }
        else
        {
          _pRevCopyProbDistr->push_back(exp(pps[copyNbr] - stot));
        }

        // define last value in probability distribution
        lastValue = (*_pRevCopyProbDistr)[copyNbr];

        // delete intermediate variables
        delete[]pps;
        pps = NULL;

        // update prior by adding next value to prior array
        if (copyNbr < _pPriorDistr->size() - 1)
        {
          (*_pPriorDistr)[copyNbr + 1] = (*_pPriorDistr)[copyNbr] / 4;
        }
        else
        {
          _pPriorDistr->push_back(((*_pPriorDistr)[copyNbr] / 4));
        }
				
				// normalisation
        stot = 0;
        for (uint j = 0; j < _pPriorDistr->size(); j++)
          stot += (*_pPriorDistr)[j];

        for (uint j = 0; j < _pPriorDistr->size(); j++)
          (*_pPriorDistr)[j] = (*_pPriorDistr)[j] / stot;

      }
      else
      {
        _pRevCopyProbDistr->push_back(0);
        lastValue = 0;
      }

      // move to next
      copyNbr++;
    }

    if (pRevPs != NULL)
      delete pRevPs;
    pRevPs = NULL;
  }

  // cleanup some mess
  delete pPs;
  pPs = NULL;

  return;
}
*/
/******************************************************************************
  Method:       UpdateCopyProbability
  Description:  adjusted to use local _priorDistrs instead of recomputing 
                all the time, the local is linked to main and can have different
                prior-input formats now.
  Author-Date:  MC-2009/10/19
  Author-Date:  MC-2013/03/11 : apply _plogNormFactors here
******************************************************************************/
void
SequenceComputation::UpdateCopyProbability(strand_modes STRAND)
{
  // if case f (fixed) : set _p(Rev)CopyProbDistr to vector[0]
  if (_priorDistrs == NULL)
  {
    // this has already been done in LinkNbrInstPrior,
    // only needs to be done once as case f will always stay fixed
    return;

  }
	
  // other cases, compute _p(Rev)CopyProbDistr based on prior.

  // reset copy probability arrays
  vector < double >::iterator iter;
  for (iter = _pCopyProbDistr->begin(); iter != _pCopyProbDistr->end(); ++iter)
  {
    (*iter) = 0;
  }
  for (iter = _pRevCopyProbDistr->begin(); iter != _pRevCopyProbDistr->end(); ++iter)
  {
    (*iter) = 0;
  }

/*
  // reset prior distribution
  for (iter = _pPriorDistr->begin(); iter != _pPriorDistr->end(); ++iter)
  {
    (*iter) = 0;
  }
	
  // initialize with 1-prior and prior
  (*_pPriorDistr)[0] = 1 - prior;
  (*_pPriorDistr)[1] = prior;
*/

	// some variables
  double nextValue, lastValue = 1;
  uint copyNbr = 1;
  vector < double >*pPs = new vector < double >;
  pPs->push_back(_dLogP0);
  ScoreVector * prior = NULL; //

  // positive strand
  if (STRAND == both || STRAND == plus_strand)
  {
    while (_length > (int) copyNbr * _wlength 
		   && lastValue >= 0.001 // leave this constraint also for bMax-case...
		   && (int)copyNbr <= _maxInst) // 
    {
      // link to the correct distribution from _priorDistrs
      prior = (*_priorDistrs)[copyNbr]->GetVector(); //
		
      // compute next value of CopyProbability
      nextValue =
        _dLogP0 - (*_plogNormFactors)[copyNbr] + ComputeNInstancesProbability(_pExpWx, copyNbr, _length,
                                               _wlength);

      //if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      if (finite(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      {
        // add element to pPs
        pPs->push_back(nextValue);

        // intermediate values
        double *pps = new double[copyNbr + 1];
        double stot = 0;
        for (uint j = 0; j < (copyNbr + 1); j++)
        {
          //pps[j] = (*pPs)[j] + log((*_pPriorDistr)[j]);
          pps[j] = (*pPs)[j] + log((*prior)[j]); //
          stot += exp(pps[j] - _dLogP0);
        }
        stot = _dLogP0 + log(stot);

        // update copy probability distribution
        for (uint j = 0; j < copyNbr; j++)
        {
          (*_pCopyProbDistr)[j] = exp(pps[j] - stot);
        }
        // add new entry to copy probability distribution
        if (copyNbr < _pCopyProbDistr->size() - 1)
        {
          (*_pCopyProbDistr)[copyNbr] = exp(pps[copyNbr] - stot);
        }
        else
        {
          _pCopyProbDistr->push_back(exp(pps[copyNbr] - stot));
        }

        // define last value in probability distribution
        lastValue = (*_pCopyProbDistr)[copyNbr];

        // delete intermediate variables
        delete[] pps;
        pps = NULL;

        /*// update prior by adding next value to prior array
        if (copyNbr < _pPriorDistr->size() - 1)
        {
          (*_pPriorDistr)[copyNbr + 1] = (*_pPriorDistr)[copyNbr] / 4;
        }
        else
        {
          _pPriorDistr->push_back(((*_pPriorDistr)[copyNbr] / 4));
        }
				
				// normalisation
        stot = 0;
        for (uint j = 0; j < _pPriorDistr->size(); j++)
          stot += (*_pPriorDistr)[j];

        for (uint j = 0; j < _pPriorDistr->size(); j++)
          (*_pPriorDistr)[j] = (*_pPriorDistr)[j] / stot;
		  */

      }
      else
      {
        _pCopyProbDistr->push_back(0);
        lastValue = 0;
      }

      // move to next
      copyNbr++;
    }
  }
	
  // cleanup
  prior = NULL; // do not delete !  

  // negative strand
  if (STRAND == both || STRAND == minus_strand)
  {
    copyNbr = 1;
    lastValue = 1;
    vector < double >*pRevPs = new vector < double >;
    pRevPs->push_back(_dLogP0);

    while (_length > (int) copyNbr * _wlength && lastValue >= 0.001
		   && (int)copyNbr <= _maxInst) // 
    {
      // link to the correct distribution from _priorDistrs
      prior = (*_priorDistrs)[copyNbr]->GetVector(); //
		
      // compute next value of CopyProbability
      nextValue =
        _dLogP0 - (*_plogNormFactors)[copyNbr] + ComputeNInstancesProbability(_pRevExpWx, copyNbr, _length,
                                               _wlength);

      // if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      if (finite(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      {
        // add element to pRevPs
        pRevPs->push_back(nextValue);

        // intermediate values
        double *pps = new double[copyNbr + 1];
        double stot = 0;
        for (uint j = 0; j < (copyNbr + 1); j++)
        {
          //pps[j] = (*pRevPs)[j] + log((*_pPriorDistr)[j]);
          pps[j] = (*pRevPs)[j] + log((*prior)[j]); //
          stot += exp(pps[j] - _dLogP0);
        }
        stot = _dLogP0 + log(stot);

        // update copy probability distribution
        for (uint j = 0; j < copyNbr; j++)
        {
          (*_pRevCopyProbDistr)[j] = exp(pps[j] - stot);
        }
        // add new entry to copy probability distribution
        if (copyNbr < _pRevCopyProbDistr->size() - 1)
        {
          (*_pRevCopyProbDistr)[copyNbr] = exp(pps[copyNbr] - stot);
        }
        else
        {
          _pRevCopyProbDistr->push_back(exp(pps[copyNbr] - stot));
        }

        // define last value in probability distribution
        lastValue = (*_pRevCopyProbDistr)[copyNbr];

        // delete intermediate variables
        delete[]pps;
        pps = NULL;

        /*// update prior by adding next value to prior array
        if (copyNbr < _pPriorDistr->size() - 1)
        {
          (*_pPriorDistr)[copyNbr + 1] = (*_pPriorDistr)[copyNbr] / 4;
        }
        else
        {
          _pPriorDistr->push_back(((*_pPriorDistr)[copyNbr] / 4));
        }
				
				// normalisation
        stot = 0;
        for (uint j = 0; j < _pPriorDistr->size(); j++)
          stot += (*_pPriorDistr)[j];

        for (uint j = 0; j < _pPriorDistr->size(); j++)
          (*_pPriorDistr)[j] = (*_pPriorDistr)[j] / stot;*/

      }
      else
      {
        _pRevCopyProbDistr->push_back(0);
        lastValue = 0;
      }

      // move to next
      copyNbr++;
    }

    if (pRevPs != NULL)
      delete pRevPs;
    pRevPs = NULL;
  }

  // cleanup some mess
  delete pPs;
  pPs = NULL;
  prior = NULL;

  return;
}


/******************************************************************************
  Method:       UpdateFixedSizeCopyProbability
  Class:        SequenceComputation
  Arguments:    int nbr, double prior, strand_modes STRAND
  
  Description:  update the probability distribution to estimate
                the number of motif instances for a fixed number of instances
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
/*void
SequenceComputation::UpdateFixedSizeCopyProbability(int nbr, double prior,
                                                    strand_modes STRAND)
{
  if (nbr < 0)
  {
    cerr <<
      "--Error-- SequenceComputation::UpdateFixedSizeCopyProbability(): number of maximal number of instances is smaller than 1."
      << endl;
    // reset local variables
    _pPriorDistr->resize(2, 0);
    _pCopyProbDistr->resize(2, 0);
    _pRevCopyProbDistr->resize(2, 0);
    return;
  }

  // define prior probability (fixed size)
  if (_pPriorDistr == NULL)
    _pPriorDistr = new vector < double >(nbr + 1, 0);

  if ((int) _pPriorDistr->size() != nbr + 1)
    _pPriorDistr->resize(nbr + 1, 0);

  // initialize with 1-prior and prior
  (*_pPriorDistr)[0] = 1 - prior;
  (*_pPriorDistr)[1] = prior;

  if (nbr > 1)
  {
    // fill in values
    double sum = 0;
    for (int i = 2; i <= nbr; i++)
    {
      (*_pPriorDistr)[i] = (*_pPriorDistr)[i - 1] / 4;
      sum += (*_pPriorDistr)[i];
    }

    // normalize
    for (int i = 0; i <= nbr; i++)
      (*_pPriorDistr)[i] = (*_pPriorDistr)[i] / sum;

  }

  // local variables
  double nextValue = 1;
  int copyNbr = 1;
  vector < double >::iterator iter;

  // positive strand
  if (STRAND == both || STRAND == plus_strand)
  {
    // initialize distribtuion
    iter = _pCopyProbDistr->begin();
    (*iter) = 1;
    iter++;
    while (iter != _pCopyProbDistr->end())
    {
      (*iter) = 0;
      iter++;
    }

    copyNbr = 1;
    double *pps = new double[nbr + 1];
    vector < double >*pPs = new vector < double >(nbr + 1, 0);
    (*pPs)[0] = _dLogP0;

    while (_length > (int) copyNbr * _wlength && copyNbr <= nbr)
    {
      // compute next value of CopyProbability
      nextValue =
        _dLogP0 + ComputeNInstancesProbability(_pExpWx, copyNbr, _length,
                                               _wlength);

      // if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      if (finite(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      {
        // add element to pPs
        (*pPs)[copyNbr] = nextValue;

        // intermediate values
        double stot = 0;
        for (int j = 0; j <= copyNbr; j++)
        {
          pps[j] = (*pPs)[j] + log((*_pPriorDistr)[j]);
          stot += exp(pps[j] - _dLogP0);
        }
        stot = _dLogP0 + log(stot);

        // update copy probability distribution
        for (int j = 0; j <= copyNbr; j++)
          (*_pCopyProbDistr)[j] = exp(pps[j] - stot);

      }
      else
      {
        (*_pCopyProbDistr)[copyNbr] = 0;
      }

      // move to next
      copyNbr++;
    }

    // delete local variables
    delete pPs;
    pPs = NULL;

    delete[]pps;
    pps = NULL;
  }

  // negative strand
  if (STRAND == both || STRAND == minus_strand)
  {
    copyNbr = 1;
    double *pps = new double[nbr + 1];
    vector < double >*pRevPs = new vector < double >(nbr + 1);
    (*pRevPs)[0] = _dLogP0;

    // initialize distribtuion
    iter = _pRevCopyProbDistr->begin();
    (*iter) = 1;
    iter++;
    while (iter != _pRevCopyProbDistr->end())
    {
      (*iter) = 0;
      iter++;
    }


    while (_length > (int) copyNbr * _wlength && copyNbr <= nbr)
    {
      // compute next value of CopyProbability
      nextValue =
        _dLogP0 + ComputeNInstancesProbability(_pRevExpWx, copyNbr, _length,
                                               _wlength);

      // if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      if (finite(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
      {
        // add element to pRevPs
        (*pRevPs)[copyNbr] = nextValue;

        // intermediate values
        double stot = 0;
        for (int j = 0; j <= copyNbr; j++)
        {
          pps[j] = (*pRevPs)[j] + log((*_pPriorDistr)[j]);
          stot += exp(pps[j] - _dLogP0);
        }
        stot = _dLogP0 + log(stot);

        // update copy probability distribution
        for (int j = 0; j <= copyNbr; j++)
          (*_pRevCopyProbDistr)[j] = exp(pps[j] - stot);

      }
      else
      {
        (*_pRevCopyProbDistr)[copyNbr] = 0;
      }

      // move to next
      copyNbr++;
    }

    if (pRevPs != NULL)
      delete pRevPs;
    pRevPs = NULL;

    // delete intermediate variables
    delete[]pps;
    pps = NULL;

  }

  return;
}
*/

/******************************************************************************
  Method:       FixCopyProbability
  Class:        SequenceComputation
  Arguments:    ScoreVector *pCopyProbValues, strand_modes STRAN
  
  Description:  Set the probability distribution to estimate the number
                of motif instances to a predefined set of values
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
/*void
SequenceComputation::FixCopyProbability(ScoreVector *pCopyProbValues, strand_modes STRAND)
{
  uint l = pCopyProbValues->size();
  if ( l != _pCopyProbDistr->size() )
  {
    _pCopyProbDistr->resize(l,0);
    _pRevCopyProbDistr->resize(l,0);
  }
  
  if ( STRAND == BOTH || STRAND == plus_strand )
  {
    for ( uint i=0; i<l; i++ )
      (*_pCopyProbDistr)[i] = (*pCopyProbValues)[i];
  }
  else if ( STRAND == minus_strand )
  {
    for ( uint i=0; i<l; i++ )
      (*_pRevCopyProbDistr)[i] = (*pCopyProbValues)[i];
  }    
  return;
}
*/

/******************************************************************************
  Method:       GetEstimatedNumberInstances
  Class:        SequenceComputation
  Arguments:    strand_modes STRAND
  
  Description:  compute the expected number of motif instances based on 
								the internal distribution _pCopyProbDistr
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
/*int
SequenceComputation::GetEstimatedNumberInstances(strand_modes STRAND)
{
  double copies = 0;
  int E = 0;

  if (STRAND == plus_strand)
  {
    for (uint j = 1; j < _pCopyProbDistr->size(); j++)
    {
      copies += j * (*_pCopyProbDistr)[j];
    }
    E = int (copies + 0.5);
  }
  else if (STRAND == minus_strand)
  {
    for (uint j = 1; j < _pRevCopyProbDistr->size(); j++)
    {
      copies += j * (*_pRevCopyProbDistr)[j];
    }
    E = int (copies + 0.5);
  }
  return E;
}
*/
/******************************************************************************
  Method:       GetNumberInstances
  Class:        SequenceComputation
  Arguments:    strand_modes STRAND
  
  Description:  compute the expected number of motif instances based on 
                the internal distribution _pCopyProbDistr
				the boolean indicates for sampling or estimating. 
  Date:         2009/10/19
  Author:       Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
int
SequenceComputation::GetNumberInstances(strand_modes STRAND)
{
  //cerr << "debug--SequenceComputation::GetNumberInstances - begin " << endl; 
  int E = 0;

  if (!_sampling) // estimation of number of instances (same as before)
  {
    double copies = 0;
    if (STRAND == plus_strand)
    {
      for (uint j = 1; j < _pCopyProbDistr->size(); j++)
        copies += j * (*_pCopyProbDistr)[j];
      E = int (copies + 0.5);
    }
    else if (STRAND == minus_strand)
    {
      for (uint j = 1; j < _pRevCopyProbDistr->size(); j++)
        copies += j * (*_pRevCopyProbDistr)[j];
      E = int (copies + 0.5);
    }
  }
  else // sampling of number of instances from _p(Rev)CopyProbDistr
  {
    // create distribution 
    Distribution *pDist = NULL;
    if (STRAND == plus_strand)
      pDist = new Distribution(_pCopyProbDistr);
    else
      pDist = new Distribution(_pRevCopyProbDistr);
    if (pDist != NULL && pDist->IsNormalized() )
       E = pDist->TakeSample();
    // cleanup
    if (pDist != NULL)
      delete pDist; 
    pDist = NULL;
  }
  //cerr << "debug--SequenceComputation::GetNumberInstances - end " << endl;
  return E;
}



/******************************************************************************
  Method:       SampleInstanceStart
  Class:        SequenceComputation
  Arguments:    vector<int> & pAlignmentVector, int n, strand_modes STRAND
  
  Description:  sample n instances from the score distribution _pExpWx. 
                The results are stored in pAlignmentVector
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::SampleInstanceStart(vector < int >&pAlignmentVector,
                                         int n, strand_modes STRAND)
{
  int ndx = -1;

  // check size of vector
  if ((int) pAlignmentVector.size() < n)
    pAlignmentVector.resize(n, 0);

  // set all values to -1 
  for (uint i = 0; i < pAlignmentVector.size(); i++)
    pAlignmentVector[i] = -1;

  // create distribution from 
  Distribution *pDist = NULL;
  if (STRAND == plus_strand)
  {
    pDist = new Distribution(_pExpWx);
  }
  else
  {
    pDist = new Distribution(_pRevExpWx);
  }

  if (pDist != NULL && pDist->IsNormalized() )
  {
    for (int i = 0; i < n; i++)
    {
      ndx = pDist->TakeSample();
      pAlignmentVector[i] = ndx;
      
      if ( n > 1 )
      {
        pDist->ApplyMask(ndx - _wlength + 1, 2 * _wlength);
        if ( !pDist->IsNormalized() )
        {
          //cerr << "SampleInstanceStart : distribution not normalized!!" << endl;
          break; 
        }
      }
      // cerr << "DEBUG Update aligment vector: " << ndx << endl;
    }

    delete pDist;
    pDist = NULL;
  }
  else
  {
    //cerr << "-- Warning -- SequenceComputation::SampleInstanceStart(): distribution is NULL." << endl;
  }
  return;
}


/******************************************************************************
  Method:       SampleUniformInstanceStart
  Class:        SequenceComputation
  Arguments:    vector<int> & pAlignmentVector, int n, strand_modes STRAND
  
  Description:  sample n instances from a uniform distribution. 
                The results are stored in pAlignmentVector
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::SampleUniformInstanceStart(vector<int>& pAlignmentVector, int n,
                                                strand_modes STRAND)
{
  int ndx = -1;

  // check size of vector
  // cerr << "Aligment vector size = " << pAlignmentVector.size() << endl;
  if ((int) pAlignmentVector.size() < n)
    pAlignmentVector.resize(n, 0);

  // set all values to -1 
  // cerr << "Initialize alignment vector" << endl;
  for (uint i = 0; i < pAlignmentVector.size(); i++)
    pAlignmentVector[i] = -1;

  // create distribution from mask vector
  // cerr << "New distribution object" <<  endl;
  Distribution *pDist = NULL;
  if (STRAND == plus_strand)
  {
    pDist = new Distribution(_pMask->GetMask());
  }
  else
  {
    pDist = new Distribution(_pRevMask->GetMask());
  }
  
  if ( pDist != NULL && pDist->IsNormalized() )
  {
    for (int i = 0; i < n; i++)
    {
      
      ndx = pDist->TakeSample();
      pAlignmentVector[i] = ndx;
      
      if ( n > 1 )
      {
        pDist->ApplyMask(ndx - _wlength + 1, 2 * _wlength);
        if ( !pDist->IsNormalized() )
        {
          //cerr << "SampleUniformInstanceStart : distribution not normalized!!" << endl;
          break;
        }
      }        
      // cerr << "DEBUG Update aligment vector: " << ndx << endl;
    }
    delete pDist;
    pDist = NULL;
  }
  else
  {
    //cerr << "-- Warning -- SequenceComputation::SampleUniformInstanceStart(): distribution is NULL." << endl;
  }
  
  return;
}


/******************************************************************************
  Method:       SampleBestInstanceStart
  Class:        SequenceComputation
  Arguments:    vector<int> & pAlignmentVector, int n, strand_modes STRAND
  
  Description:  Select the n best instances from the score distribution _pExpWx. 
                The results are stored in pAlignmentVector
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::SelectBestInstanceStart(vector < int >&pAlignmentVector,
                                             int n, strand_modes STRAND)
{
  int ndx = -1;

  // check size of vector
  if ((int) pAlignmentVector.size() < n)
    pAlignmentVector.resize(n, 0);

  // set all values to -1 
  for (uint i = 0; i < pAlignmentVector.size(); i++)
    pAlignmentVector[i] = -1;

  // create distribution from 
  Distribution *pDist = NULL;
  if (STRAND == plus_strand)
  {
    pDist = new Distribution(_pExpWx);
  }
  else
  {
    pDist = new Distribution(_pRevExpWx);
  }

  if ( pDist != NULL  && pDist->IsNormalized() )
  {
    for (int i = 0; i < n; i++)
    {
      ndx = pDist->SelectMax();
      pAlignmentVector[i] = ndx;
      
      if ( n > 1 )
      {
        pDist->ApplyMask(ndx - _wlength + 1, 2 * _wlength);
        if ( !pDist->IsNormalized() )
        {
          //cerr << "SelectBestInstanceStart : distribution not normalized!!" << endl;
          break;
        }
      }        
    }

		if ( pDist != NULL )
			delete pDist;
    pDist = NULL;
  }
  else
  {
    //cerr << "-- Warning -- SequenceComputation::SelectBestInstanceStart(): distribution is NULL." << endl;
  }
  return;
}



/******************************************************************************
  Method:       GetWxAt
  Class:        SequenceComputation
  Arguments:    int ndx, strand_modes STRAND
  
  Description:  get the value Wx of the instance at position ndx
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
SequenceComputation::GetWxAt(int ndx, strand_modes STRAND)
const
{
	if (ndx < 0 || ndx >= _length)
  {
    return 0;
  }
  else if (STRAND == plus_strand)
  {
    return (*_pExpWx)[ndx];
  }
  else if (STRAND == minus_strand)
  {
    return (*_pRevExpWx)[ndx];
  }
	else
	{  
		return 0;
	}
}

/******************************************************************************
  Method:       GetMotifScoreAt
  Class:        SequenceComputation
  Arguments:    int ndx, strand_modes STRAND
  
  Description:  get the motif score of the instance at position ndx
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
SequenceComputation::GetMotifScoreAt(int ndx, strand_modes STRAND)
const
{
  if (ndx < 0 || ndx >= _length)
  {
    return 0;
  }
  else if (STRAND == plus_strand)
  {
    return (*_pMotifScore)[ndx];
  }
  else if (STRAND == minus_strand)
  {
    return (*_pRevMotifScore)[ndx];
  }
	else
	{
		return 0;
	}
}
/******************************************************************************
  Method:       GetPspScoreAt
  Description:  get the PSP score of the instance at position ndx
  
  Date:         2012/09/19
  Author:       Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
double
SequenceComputation::GetPspScoreAt(int ndx, strand_modes STRAND)
const
{
  if (ndx < 0 || ndx >= _length)
  {
    return 0;
  }
  else if (STRAND == plus_strand)
  {
    return (*_pPspScore)[ndx];
  }
  else if (STRAND == minus_strand)
  {
    return (*_pRevPspScore)[ndx];
  }
	else
	{
		return 0;
	}
}
/******************************************************************************
  Method:       GetBackgroundScoreAt
  Class:        SequenceComputation
  Arguments:    int ndx, strand_modes STRAND
  
  Description:  Get the background score of the instance at position ndx
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
SequenceComputation::GetBackgroundScoreAt(int ndx, strand_modes STRAND)
const
{
  if (ndx < 0 || ndx >= _length)
  {
    return 0;
  }
  else if (STRAND == plus_strand)
  {
    return (*_pBackgroundScore)[ndx];
  }
  else if (STRAND == minus_strand)
  {
    return (*_pRevBackgroundScore)[ndx];
  }
	else
	{
		return 0;
	}
}

/******************************************************************************
  Method:       GetLogBackgroundScore
  Class:        SequenceComputation
  Arguments:    strand_modes STRAND
  
  Description:  get the log background score of the complete sequence
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
SequenceComputation::GetLogBackgroundScore(strand_modes STRAND)
const
{
  if (STRAND == plus_strand)
  {
    return _dLogP0;
  }
  else if (STRAND == minus_strand)
  {
    return _dRevLogP0;
  }
  else
  {
    return 0;
  }
}

/******************************************************************************
  Method:       UpdateMask
  Class:        SequenceComputation 
  Arguments:    int start, int w, int value, strand_modes STRAND
  
  Description:  set the mask to value at the positions start->start+w-1
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::UpdateMask(int start, int w, int value,
                                strand_modes STRAND)
{
  if (STRAND == plus_strand)
  {
    _pMask->UpdateMask(start, w, value);
  }
  else if (STRAND == minus_strand)
  {
    _pRevMask->UpdateMask(start, w, value);
  }
  return;
}


/******************************************************************************
  Method:       GetMask
  Class:        SequenceComputation
  Arguments:    strand_modes STRAND
  
  Description:  get a pointer to the internal masking vector
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
MaskVector *
SequenceComputation::GetMask(strand_modes STRAND)
{
  if (STRAND == plus_strand)
  {
    return _pMask;
  }
  else if (STRAND == minus_strand)
  {
    return _pRevMask;
  }
  else
  {
    return NULL;
  }
}


/******************************************************************************
  Method:       GetCopyProbabilityAt
  Class:        SequenceComputation
  Arguments:    int nbr, strand_modes STRAND
  
  Description:  get the probability of finding nbr instances on STRAND
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
SequenceComputation::GetCopyProbabilityAt(int nbr, strand_modes STRAND)
{
  if (STRAND == plus_strand && nbr >= 0
      && nbr <= (int) _pCopyProbDistr->size())
  {
    return (*_pCopyProbDistr)[nbr];
  }
  else if (STRAND == minus_strand && nbr >= 0
           && nbr <= (int) _pRevCopyProbDistr->size())
  {
    return (*_pRevCopyProbDistr)[nbr];
  }
  else
  {
    return 0;
  }
}

/******************************************************************************
  Method:       SetMotifLength
  Class:        SequenceComputation
  Arguments:    int wLength
  
  Description:  set the motif length and update the internal mask 
                !!---
                  if you use this function, make sure that you also update
								  the background scores since they depend on the motif length
                ---!!
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::SetMotifLength(int wLength)
{
  if ( _wlength != wLength)
  {
    // cerr << "DEBUG define new motif length" << wLength << " from " << _length << endl;
    _wlength = wLength;

    // update mask (last wLength - 1 positions set to zero);
    _pMask->UpdateMask(_length - _wlength + 1, _wlength, 0);
    _pRevMask->UpdateMask(_length - _wlength + 1, _wlength, 0);
    int start = _length - _wlength +1;
    if (start < 0) {start = 0;}
    for (int i = start; i < _length; i++)
    {
      (*_pPspScore)[i] = 0;
      (*_pRevPspScore)[i] = 0;
    }
  }
  return;
}


/******************************************************************************
  Method:       LogLikelihoodScore
  Class:        SequenceComputation
  Arguments:    vector<int> * pAlign
  
  Description:  compute the loglikelihood score given the alignment positions
  
  Date:         2003/07/18
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double 
SequenceComputation::LogLikelihoodScore(vector<int> * pAlign, strand_modes STRAND)
{
  double score = 0, cp = 0;
  int i = 0;
  int nbr = _pCopyProbDistr->size();

  if (STRAND == plus_strand)
  {
    score = (*_pCopyProbDistr)[0];
    for (i = 0; i < (int) pAlign->size(); i++)
    {
      cp = 0;
      if (i + 1 < nbr)
        cp = SumArray(_pCopyProbDistr, i + 1, nbr - i - 1);
      cp = cp * (*_pExpWx)[(*pAlign)[i]];
      score += cp;
    }
  }
  else if (STRAND == minus_strand)
  {
    score = (*_pRevCopyProbDistr)[0];
    for (i = 0; i < (int) pAlign->size(); i++)
    {
      cp = 0;
      if (i + 1 < nbr)
        cp = SumArray(_pRevCopyProbDistr, i + 1, nbr - i - 1);
      cp = cp * (*_pRevExpWx)[(*pAlign)[i]];
      score += cp;
    }
  }

  return log(score);
}

/******************************************************************************
  Description:  set the PSP score of all segments
                and update the number of instances prior distribution !! 
	            and also compute the prior-NORM-factor, store in _pNormFactors
  Date:         2013/03/10
  Author:       Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
void
SequenceComputation::UpdatePspScores(ScoreVector * psp, strand_modes STRAND)
{
  // define some iterators
  vector < double >::iterator iter1;
  int L = psp->size();
  if (L != _length + 1) // should normally be ok (checked before)
  {
    cerr << "--Error--SeqComp::UpdatePspScores(): inconsistent psp/seqlengths (" 
		  << L << "," << _length << ")." << endl; 
    //return; // no return, it will result in an error to be fixed
  }
  // first update the prior on the NUMBER of instances in this sequence
  // i.e. Pr(c=0) == PSP[0] and Pr(c=1) == 1-PSP[0]
    // for c > 1 :recaculate proportional probabilities. 
  // this means overwriting on the earlier entries of vector<Distribution*> * _priorDistrs; 
  // so you can have different Number-priors for different sequences now !
  // REPORT changes (compared to -p input) to the user
  // if -p was a fixed input (eg f0_1), then there is NO update 
  // -> this is the case if _priorDistrs == NULL => then do nothing
  double prior0 = (*psp)[0];
  double prior1 = 1- prior0;
  vector<double> * distr = NULL;
  double norm, sum;
  if (_priorDistrs != NULL) 
  {
    cerr << "WARNING: overwrite prior distribution on NUMBER of instances for sequence [" 
         << *_pSeqObj->GetID() << "]: "<< endl;
    // update each distribution ( c= 1, ... c= M)
    for (int c = 1; c < int(_priorDistrs->size()); c++)
    {
      distr = (*_priorDistrs)[c]->GetVector();
      cerr << "[c=" << c << "]::Before(";
      for(int cc = 0; cc < int(distr->size()); cc++)
		   { cerr << (*distr)[cc] << " ";}
      norm = (*distr)[0]  + (*distr)[1];
      // overwrite prior0 and prior1
      (*distr)[0] = prior0; (*distr)[1] = prior1;
      // normalize and assign higher probs
      sum = norm;
      for (int cc = 2; cc < int(distr->size()); cc++)
      { (*distr)[cc] /= norm; sum += (*distr)[cc];}
      for (int cc = 0; cc < int(distr->size()); cc++)
      { (*distr)[cc] /= sum;}
      cerr << "; After(";
      for (int cc = 0; cc < int(distr->size()); cc++)
      { cerr << (*distr)[cc] << " ";}
      cerr << endl;
    }
  }

  // update each PSP prior segment score 
  // and renormalize so the sum equals 1 
  // plus strand
  if (STRAND == both || STRAND == plus_strand)
  { sum = 0;
    for (int i = 0; i < _length; i++)
    {
      if (i > _length - _wlength) {(*_pPspScore)[i] = 0;}
      else {(*_pPspScore)[i] = (*psp)[i+1];}
      sum += (*_pPspScore)[i];
    }
    for (int i = 0; i < _length; i++)
    { (*_pPspScore)[i] /= sum;}
  }
  // minus strand // apply same probs on the reverse strand 
  if (STRAND == both || STRAND == minus_strand)
  { sum = 0;
    for (int i = 0; i < _length; i++)
    {
      if (i > _length - _wlength) {(*_pRevPspScore)[i] = 0;}
      else {(*_pRevPspScore)[i] = (*psp)[L-i-1];}
      sum += (*_pRevPspScore)[i];
    }
    for (int i = 0; i < _length; i++)
    { (*_pRevPspScore)[i] /= sum;}
  }

  // now also compute/store the NORM-factor for c=1->M
  // NORM-factor = sum of all possible combined site-priors of sets of sites
  // (code based on utilities->ComputeNInstancesProbability - omit "C" in there
  for(int c = 1; c < _maxInst+1; c++)
  {
    (*_plogNormFactors)[c] = INCLUSIVE::ComputeNInstancesProbability(
        _pPspScore, c, _length, _wlength);
  }
  // (*_plogNormFactors)[1] should be zero by this computation (=log(1));
	  cerr << "debug-followup-SequenceComputation::UpdatePSPScores:" << 
		  "_plogNormFactors psp =(";
	  for (int i = 0; i < _maxInst+1; i++)
		  cerr << (*_plogNormFactors)[i] << " ";
	  cerr << ")" << endl;
  return;
}

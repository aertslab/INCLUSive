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
  _pPrior2NormFactor = new ScoreVector(10,0);
  _pRevPrior2NormFactor = new ScoreVector(10,0);
  //_pPriorDistr = new ScoreVector(10, 0);
  _priorDistrs = NULL; // further in LinkNbrInstPrior

  // initial motif length is zero
  _wlength = 0;

  // initially the format is correct
  _correct = true;
  PSP_s = false; //  default no PSP sampling impact

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
	
  if (_pPrior2NormFactor != NULL)
    delete _pPrior2NormFactor;
  _pPrior2NormFactor = NULL;

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

  if (_pRevPrior2NormFactor != NULL)
    delete _pRevPrior2NormFactor;
  _pRevPrior2NormFactor = NULL;

  if (_pRevExpWx != NULL)
    delete _pRevExpWx;
  _pRevExpWx = NULL;

  if (_pRevMask != NULL)
    delete _pRevMask;
  _pRevMask = 0;

  if (_pRevPspScore != NULL)
    delete _pRevPspScore;
  _pRevPspScore = NULL;

  _priorDistrs = NULL; // do not delete, is done in main

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
  Description:  LoadNbrInstFixed : copy values into_p(Rev)CopyProb fixed
  Author_Date:  MC_2013/12/30
******************************************************************************/
void 
SequenceComputation::LoadNbrInstFixed(ScoreVector * probvector)
{ 
  //cerr << "debug--SequenceComputation::LoadNbrInstFixed - begin" << endl;
  // set the maximum searchspace
  _maxInst = (int)probvector->size()-1;
  if (_maxInst > _length/_wlength ) // reset max for short sequences
    _maxInst = (int)(_length/_wlength);

  _priorDistrs = NULL; // not used in case f

  // copy the entries to the fixed copy-distribution
  for(int i = 0; i < (int)probvector->size(); i++)
  {
    if (i <= 9) // replace (see constructor new ScoreVector(10, 0);)
    {
      (*_pCopyProbDistr)[i] = (*probvector)[i];
      (*_pRevCopyProbDistr)[i] = (*probvector)[i];
    }
    else // add
    {
      _pCopyProbDistr->push_back((*probvector)[i]);
      _pRevCopyProbDistr->push_back((*probvector)[i]);
    }
  }
  //cerr << "debug--SequenceComputation::LoadNbrInstFixed - end" << endl;
  return;
}
/******************************************************************************
  Description:  LoadNbrInstPrior : link pointer and set _maxInst
  Author_Date:  MC_2013/12/30
******************************************************************************/
void 
SequenceComputation::LinkNbrInstPrior(vector<Distribution*>* distrs)
{ 
  // set the maximum searchspace
  _maxInst = (int)distrs->size();
  if (_maxInst > _length/(_wlength)) // reset max for short sequences
    _maxInst = (int)(_length/(_wlength));

  // link
  _priorDistrs = distrs; 

  //cerr << "debug--SequenceComputation::LinkNbrInstPrior - end" << endl;
  return;
}

/******************************************************************************
  Method:       UpdateInstanceExpWx
  Class:        SequenceComputation
  Arguments:    strand_modes STRAND
  
  Description:  set the score of all instances based on the resp. motif scores
                and the background scores
                // the scores are absolute (no logarithms)  
  Date:         2012/09/19
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
	            Revised: Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
void
SequenceComputation::UpdateInstanceExpWx(strand_modes STRAND)
{
  //cerr << "SeqComp::UpdateInstanceExpWx - begin" << endl;
  // define some iterators
  vector < double >::iterator iter1;
  vector < double >::iterator iter2;
  vector < double >::iterator iter3;
  vector < int >::iterator iter4;
  // plus strand
  if (STRAND == both || STRAND == plus_strand)
  {
    iter1 = _pMotifScore->begin();
    iter2 = _pBackgroundScore->begin();
    iter3 = _pExpWx->begin();
    iter4 = _pMask->begin();

    while (iter1 != _pMotifScore->end())
    {
      (*iter3) = exp((*iter1) - (*iter2)) * (*iter4);// no LOG !!
      //cerr << (*iter3) << " ";
      // augment iterators
      iter1++;
      iter2++;
      iter3++;
      iter4++;
    }
  }
  // minus strand
  if (STRAND == both || STRAND == minus_strand)
  {
    iter1 = _pRevMotifScore->begin();
    iter2 = _pRevBackgroundScore->begin();
    iter3 = _pRevExpWx->begin();
    iter4 = _pRevMask->begin();
    while (iter1 != _pRevMotifScore->end())
    {
      (*iter3) = exp((*iter1) - (*iter2)) * (*iter4);
      // augment iterators
      iter1++;
      iter2++;
      iter3++;
      iter4++; 
    }
  }
  //cerr << "SeqComp::UpdateInstanceExpWx - end" << endl;
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
  Description:  computation of Pr(Qk=c|Sk,motif,bg) for c = 0,...,Max
                adjusted to use local _priorDistrs (as prior1) instead of recomputing 
                all the time, the local is linked to main and can have different
                prior-input formats now.
  Author-Date:  MC-2009/10/19
  Author-Date:  MC-2013/12/28 : apply 'prior2' based on _p(rev)PSPScore
                and compute norm-cte by means of ComputeNInstancesProbability
                (accounts also for computation of "C" in old case)
******************************************************************************/
void
SequenceComputation::UpdateCopyProbability(strand_modes STRAND)
{
  // if case f (fixed) 
  if (_priorDistrs == NULL)
  {
    // _p(rev)CopyProbDistr has been set fixed in LoadNbrInstFixed
    return;
  }


  // other cases, compute _p(Rev)CopyProbDistr based on prior.
  // first compute prior2 normalization constantes for c=1->M
  // this is only needed once (first time this function is called)
  if ( (*_pPrior2NormFactor)[0] == 0) // first time
  { 
    for(int c = 1; c <= _maxInst; c++)
    {
      if (c <= 10)
      {
        (*_pPrior2NormFactor)[c-1] = 
          ComputeNInstancesProbability(_pPspScore,c,_length,_wlength);
        (*_pRevPrior2NormFactor)[c-1] = 
          ComputeNInstancesProbability(_pRevPspScore,c,_length,_wlength);
      }
      else
      {
        _pPrior2NormFactor->push_back(
          ComputeNInstancesProbability(_pPspScore,c,_length,_wlength));
        _pRevPrior2NormFactor->push_back( 
          ComputeNInstancesProbability(_pRevPspScore,c,_length,_wlength));
      }
    }
	  /* // report for evaluation
    for(int c = 1; c <= _maxInst; c++)
    {
      Check = log((double) (_length - (c*_wlength) +1));
      for ( int i = 2; i <= (int)c; i++)
        Check += log((double)(_length-(c*_wlength) + i))- log((double) i);
      cerr << "Temporary-Info::for c=" << c << ":old-C=" << Check << "/";
      cerr << "ComputeNInstancesProbability-C=" << (*_pPrior2NormFactor)[c-1];
      cerr << ",RevComputeNInstancesProbability-C=" << (*_pRevPrior2NormFactor)[c-1] << endl;
    }
	*/
  }

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
  // pPs vector to store log[Pr(Sk|Akc)*Pr(Akc)] (terms without prior1)
  vector < double >*pPs = new vector < double >;
  pPs->push_back(_dLogP0);
  ScoreVector * prior1 = NULL; //
  ScoreVector * pSiteScore = new ScoreVector(_length, 0);
  
  // positive strand
  if (STRAND == both || STRAND == plus_strand)
  {
    // compute the motif-psp score for each segment
    for(int p = 0; p < _length; p++)
    { (*pSiteScore)[p] = (*_pExpWx)[p] * (*_pPspScore)[p];}

    // compute for sequentially higher copyNbr until lastvalue < 0.001 or Max
    while (_length > (int) copyNbr * _wlength 
		   && lastValue >= 0.001 // leave this constraint also for bMax-case...
		   && (int)copyNbr <= _maxInst) // 
    {
      // link to the correct 'prior1' distribution from _priorDistrs
      prior1 = (*_priorDistrs)[copyNbr-1]->GetVector(); // '-1' 2013/07/23

      // compute next value of CopyProbability
      nextValue =  _dLogP0 + 
            ComputeNInstancesProbability(pSiteScore, copyNbr, _length, _wlength)
            - (*_pPrior2NormFactor)[copyNbr-1];

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
          pps[j] = (*pPs)[j] + log((*prior1)[j]); //
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
  prior1 = NULL; // do not delete !  

  // negative strand
  if (STRAND == both || STRAND == minus_strand)
  {
    // compute the motif-psp score for each segment
    for(int p = 0; p < _length; p++)
    { (*pSiteScore)[p] = (*_pRevExpWx)[p] * (*_pRevPspScore)[p];}

    copyNbr = 1;
    lastValue = 1; nextValue = 1;
    vector < double >*pRevPs = new vector < double >;
    pRevPs->push_back(_dLogP0);

    while (_length > (int) copyNbr * _wlength && lastValue >= 0.001
		   && (int)copyNbr <= _maxInst) // 
    {
      // link to the correct 'prior1' distribution from _priorDistrs
      prior1 = (*_priorDistrs)[copyNbr-1]->GetVector(); // '-1' 2013/07/23
      
      // compute next value of CopyProbability
      nextValue = _dLogP0 
            + ComputeNInstancesProbability(pSiteScore,copyNbr, _length,
                                               _wlength)
            - (*_pRevPrior2NormFactor)[copyNbr-1];
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
          pps[j] = (*pRevPs)[j] + log((*prior1)[j]); //
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
  delete pPs; pPs = NULL;
  prior1 = NULL;
  if (pSiteScore != NULL) delete pSiteScore; pSiteScore = NULL;

  return;
}

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
SequenceComputation::GetNumberInstances(strand_modes STRAND, bool bSelect)
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
    { if (bSelect) { E = pDist->SelectMax();}
      else { E = pDist->TakeSample();}
    }
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
    Date:         2013/12/30 : add impact PSP sampling
  Author:       Marleen Claeys 
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
  ScoreVector * pSiteScore = NULL;
  if (STRAND == plus_strand)
  {
    if (PSP_s)
    { // update motif-scores to motif-psp-scores (test "extra" impact "sampling")
      pSiteScore = new ScoreVector(_length, 0);
      for(int p = 0; p < _length; p++)
      { (*pSiteScore)[p] = (*_pExpWx)[p] * (*_pPspScore)[p];}
    } 
    else { pSiteScore = _pExpWx;}
  }
  else
  {
    if (PSP_s)
    { // update motif-scores to motif-psp-scores (test "extra" impact "sampling")
      pSiteScore = new ScoreVector(_length, 0);
      for(int p = 0; p < _length; p++)
      { (*pSiteScore)[p] = (*_pRevExpWx)[p] * (*_pRevPspScore)[p];}
    } 
    else { pSiteScore = _pRevExpWx;}
  }
  // create distribution
  pDist = new Distribution(pSiteScore);

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
  if (PSP_s && pSiteScore != NULL) {delete pSiteScore;} 
  pSiteScore = NULL;

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
  Method:       SelectBestInstanceStart
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
  ScoreVector * pSiteScore = NULL;
  if (STRAND == plus_strand)
  {
    if (PSP_s)
    { // update motif-scores to motif-psp-scores (test "extra" impact "sampling")
      pSiteScore = new ScoreVector(_length, 0);
      for(int p = 0; p < _length; p++)
      { (*pSiteScore)[p] = (*_pExpWx)[p] * (*_pPspScore)[p];}
    } 
    else { pSiteScore = _pExpWx;}
  }
  else
  {
    if (PSP_s)
    { // update motif-scores to motif-psp-scores (test "extra" impact "sampling")
      pSiteScore = new ScoreVector(_length, 0);
      for(int p = 0; p < _length; p++)
      { (*pSiteScore)[p] = (*_pRevExpWx)[p] * (*_pRevPspScore)[p];}
    } 
    else { pSiteScore = _pRevExpWx;}
  }
  // create distribution
  pDist = new Distribution(pSiteScore);

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
    delete pDist;
    pDist = NULL;
  }
  else
  {
    //cerr << "-- Warning -- SequenceComputation::SelectBestInstanceStart(): distribution is NULL." << endl;
  }
  if (PSP_s && pSiteScore != NULL){ delete pSiteScore;} 
  pSiteScore = NULL;

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
  Date:         2013/12/28
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
  Date:         2014/01/14 (debug first reset mask)
  Author:       Marleen Claeys
******************************************************************************/
void
SequenceComputation::SetMotifLength(int wLength)
{
  if ( _wlength != wLength)
  {
    // cerr << "DEBUG define new motif length" << wLength << " from " << _length << endl;
    _wlength = wLength;

    // update mask (last wLength - 1 positions set to zero);
    _pMask->ResetMask();
    _pMask->UpdateMask(_length - _wlength + 1, _wlength, 0);
    _pRevMask->ResetMask();
    _pRevMask->UpdateMask(_length - _wlength + 1, _wlength, 0);

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
	//cerr << "DEBUG: SeqComp::LogLikelihoodScore - begin" << endl;
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
	  last L-w values are already set to zero (in both strands if so)
	  L entries : for + strand, symmetrical if - strand needed
	  2L entries : for + respectively - strand (if needed)
	  convert scores to S(x)/(1-S(x))
	  NO NORMALIZATION HERE 
  Date:         2014/01/12
  Author:       Marleen Claeys <mclaeys@qatar.net.qa>
  
******************************************************************************/
void
SequenceComputation::LoadPspScores(ScoreVector * psp, strand_modes STRAND, bool psp_s)
{
//cerr <<"debug--SequenceComputation::LoadPspScores-begin" << endl;
  // set local _PSPimpact
  PSP_s = psp_s;
//cerr << "seq=" << *_pSeqObj->GetID() << endl;
  // define some iterators
  //vector < double >::iterator iter1;
  int L = psp->size();
  // update each PSP prior score
	// Mark : the input psp values are SCORES, not a distribution
	// Mark : the values in psp have been rescaled [0.1-0.9] whilst reading from file
	// Mark : last 'w' psp-positions were reset to zero
  // plus strand
  if (STRAND == both || STRAND == plus_strand)
  { 
    for (int i = 0; i < _length; i++)
    { (*_pPspScore)[i] = (*psp)[i]; } 
    if ((*psp)[0] != 1) // else all supplied entries were uniform
    { for (int i = 0; i <= _length - _wlength; i++)
      { // convert the S(x) score into Pr(x) probability 
        // last 'w' values are not relevant (skip)
        // Pr(x) ~ S(x)/(1-S(x))
        (*_pPspScore)[i] /= (1-(*psp)[i]); 
      }
    }
  }
  // minus strand 
  if (STRAND == both || STRAND == minus_strand)
  {
    if (L == 2*_length)
    { 
      for (int i = 0; i < _length; i++)
      { (*_pRevPspScore)[i] = (*psp)[_length+i]; }
      if ((*psp)[_length] != 1)
      { for (int i = 0; i <= _length - _wlength; i++)
        { // convert the S(x) score into Pr(x) probability 
          (*_pRevPspScore)[i] /= (1-(*psp)[_length+i]);
        }
      }
    }
    else // use "shift-symmetrical" values of + strand
    { 
      for (int i = 0; i <= _length - _wlength; i++)
      { (*_pRevPspScore)[i] = (*psp)[_length - _wlength - i];}
      for (int i = _length - _wlength +1; i < _length; i++)
      { (*_pRevPspScore)[i] = 0;}
      if ((*psp)[_length - _wlength] != 1)
      {
        for (int i = 0; i <= _length - _wlength; i++)
        { // convert the S(x) score into Pr(x) probability 
          (*_pRevPspScore)[i] /= (1-(*psp)[_length - _wlength - i]);
        }
      }
    }
  }
  //cerr <<"debug--SequenceComputation::LoadPspScores-end" << endl;
  return;
}

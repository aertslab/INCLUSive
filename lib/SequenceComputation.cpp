#include "SequenceComputation.h"
#include "utilities.h"
#include <cmath>

using namespace INCLUSIVE;

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
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

  _pRevMotifScore = new ScoreVector(_length);
  _pRevBackgroundScore = new ScoreVector(_length, 0);
  _pRevExpWx = new ScoreVector(_length, 0);
  _pRevMask = new MaskVector(_length);

  // set initial length of copy probability to 10
  _pCopyProbDistr = new ScoreVector(10, 0);
  _pRevCopyProbDistr = new ScoreVector(10, 0);
  _pPriorDistr = new ScoreVector(10, 0);

  // initial motif length is zero
  _wlength = 0;

  // initially the format is correct
  correct = true;

}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
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

  // other 
  if (_pPriorDistr != NULL)
    delete _pPriorDistr;
  _pPriorDistr = NULL;

}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::UpdateInstanceExpWx(strand_modes STRAND)
{
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
      (*iter3) = exp((*iter1) - (*iter2)) * (*iter4);
      // cerr << (*iter3) << " ";
      // augment iterators
      iter1++;
      iter2++;
      iter3++;
      iter4++;
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
  return;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::UpdateCopyProbability(double prior, strand_modes STRAND)
{
  // reset copy probability arrays
  vector < double >::iterator iter;
  for (iter = _pCopyProbDistr->begin(); iter != _pCopyProbDistr->end();
       ++iter)
  {
    (*iter) = 0;
  }
  for (iter = _pRevCopyProbDistr->begin(); iter != _pRevCopyProbDistr->end();
       ++iter)
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

      if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
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

      if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
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


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
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

      if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
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

      if (!isinf(nextValue) && !isnan(nextValue) && nextValue != -HUGE_VAL)
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


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
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

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
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



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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

  if (pDist != NULL)
  {
    for (int i = 0; i < n; i++)
    {
      ndx = pDist->TakeSample();
      pAlignmentVector[i] = ndx;
      
      if ( n > 1 )
      {
        pDist->ApplyMask(ndx - _wlength + 1, 2 * _wlength);
        if ( pDist == NULL )
          break; 
      }
      // cerr << "DEBUG Update aligment vector: " << ndx << endl;
    }

    delete pDist;
    pDist = NULL;
  }
  else
  {
    cerr << "-- Warning -- SequenceComputation::SampleInstanceStart(): distribution is NULL." << endl;
  }
  return;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
  
  if ( pDist != NULL )
  {
    for (int i = 0; i < n; i++)
    {
      ndx = pDist->TakeSample();
      pAlignmentVector[i] = ndx;
      
      if ( n > 1 )
      {
        pDist->ApplyMask(ndx - _wlength + 1, 2 * _wlength);
        if ( pDist == NULL )
          break;
      }        
      // cerr << "DEBUG Update aligment vector: " << ndx << endl;
    }
    delete pDist;
    pDist = NULL;
  }
  else
  {
    cerr << "-- Warning -- SequenceComputation::SampleUniformInstanceStart(): distribution is NULL." << endl;
  }
  
  return;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
SequenceComputation::SelectBestInstanceStart(vector < int >&pAlignmentVector,
                                             int n, strand_modes STRAND)
{
  int ndx = -1;

  // cerr << "DEGUG start selection " << endl;

  // check size of vector
  if ((int) pAlignmentVector.size() < n)
    pAlignmentVector.resize(n, 0);

  // set all values to -1 
  for (uint i = 0; i < pAlignmentVector.size(); i++)
    pAlignmentVector[i] = -1;

  // cerr << "DEGUG alignment vector initialized " << endl;

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

  if ( pDist != NULL )
  {
    for (int i = 0; i < n; i++)
    {
      ndx = pDist->SelectMax();
      pAlignmentVector[i] = ndx;
      
      if ( n > 1 )
      {
        pDist->ApplyMask(ndx - _wlength + 1, 2 * _wlength);
        if ( pDist == NULL )
          break;
      }        
      // cerr << "DEBUG Alignmentvector index " << ndx << " number " << i << endl;
    }
    delete pDist;
    pDist = NULL;
  }
  else
  {
    cerr << "-- Warning -- SequenceComputation::SelectBestInstanceStart(): distribution is NULL." << endl;
  }
  return;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
SequenceComputation::GetWxAt(int ndx, strand_modes STRAND)
{
  double wx = 0;

  if (ndx < 0 || ndx > _length)
  {
    wx = 0;
  }
  else if (STRAND == plus_strand)
  {
    wx = (*_pExpWx)[ndx];
  }
  else if (STRAND == minus_strand)
  {
    wx = (*_pRevExpWx)[ndx];
  }

  return wx;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
SequenceComputation::GetLogBackgroundScore(strand_modes STRAND)
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
  }
  return;
}


/******************************************************************************
  Method: LogLikelihoodScore(vector<int> *pAlign)
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/07/18
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
 
  SequenceComputation::LogLikelihoodScore(vector < int >*pAlign,
                                          strand_modes STRAND)
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

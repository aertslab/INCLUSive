#include "utilities.h"
#include <math.h>
namespace INCLUSIVE {

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
  int Site2Index(SequenceObject * pSeq, strand_modes s, int start, int order) {
    int index = 0;
    int nt = 0;

    if (order < 0) {
      cerr <<
        "Error: Site2index   order set to value smaller than 0, returning -1."
        << endl;
      return -1;
    } if (start < 0) {
      cerr <<
        "Error: Site2index   start set to value smaller than 0, returning -1."
        << endl;
      return -1;
    }

    for (int j = 0; j < order; j++)
    {
      nt = pSeq->GetNucleotideAt(s, start + j);
      if (nt == -1)
      {
        // this is a non ACGT symbol
        index = -1;
        break;
      }
      index += (((int) pow(4.0, order - j - 1)) * nt);
    }
    return index;
  }

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
  double SumArray(const double *pt, const int start, const int number) {
    int i = 0;
    double sum = 0;
    for (i = 0; i < number; i++)
    {
      sum += *(pt + start + i);
    }
    return sum;
  }

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
  double SumArray(vector < double >::iterator pt,
                  const int start, const int number) {
    int i = 0;
    double sum = 0;
    for (i = 0; i < number; i++)
    {
      sum += *(pt + start + i);
    }
    return sum;
  }


  double SumArray(vector < double >*pt, const int start, const int number) {
    int i = 0;
    double sum = 0;

    for (i = 0; i < number; i++)
    {
      sum += (*pt)[start + i];
    }
    return sum;
  }


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
  double LogGamma(const double x){
    static const double coefficients[6] = {
      76.18009172947146,
      -86.50532032941677,
      24.01409824083091,
      -1.231739572450155,
      0.1208650973866179e-2,
      -0.5395239384953e-5
    };
    double 
      t1 = x,
      t2 = x,
      t3 = 1.000000000190015,
      value = 0;
    
    t2 += 5.5;
    t2 -= (x+0.5)* log(t2);
    for (int j=0; j<6; j++)
      t3 += coefficients[j]/(++t1);
    
    value = -t2 + log(2.5066282746310005 * t3 / x);    
    return value;
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
    LogNormDirichlet(vector<double>* pVec)
  {
    int l = (int)pVec->size();
    double value = 0;
    double sum = SumArray(pVec,0,l);
    for (int j=0; j<l; j++)
      value += LogGamma((*pVec)[j]);
    value -= LogGamma(sum);
    
    return value;
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
    LogDirichlet(vector<double>* pParam, vector<double>* pVec)
  {
    // get normalization factor
    double normC = LogNormDirichlet(pParam);
    double value = normC;
    int l1 = (int)pVec->size();
    int l2 = (int)pParam->size();
    
    if ( l1 != l2 )
    {
      cerr << "--Error-- INCLUSive::LogDirichlet: Dimension mismatch." << endl;
      return -HUGE_VAL;
    }
    
    for (int j=0; j<l1; j++)
      value += ((*pParam)[j] - 1) * log((*pVec)[j]);

    return value;
  }
  
/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
  int MaxIndex(vector < double >*pt) {
    double m = (*pt)[0];
    int ndx = 0;
    for (uint i = 0; i < pt->size(); i++)
    {
      if ((*pt)[i] > m)
      {
        ndx = i;
        m = (*pt)[ndx];
      }
    }
    return ndx;
  }

  
/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
  void SegmentLogMatrixScore(ScoreVector & pmx,
                             SequenceObject * pSeq,
                             strand_modes s, PWM * pMatrix) {

    // check length
    int W = pMatrix->Length();
    if (pSeq->Length() < W)
    {
      //pmx = NULL;
      return;
    }
    // total number of motif positions
    int Lw = pSeq->Length() - W + 1;

    // iterator
    vector < double >::iterator fIter = pmx.begin();
    double value = 0;
    // cerr << "logPmx = ";
    for (int i = 0; i < Lw; i++)
    {
      double v = 0;
      for (int j = 0; j < W; j++)
      {
        int nt = pSeq->GetNucleotideAt(s, i + j);
        if (nt != -1)
        {
          value = pMatrix->GetValueAt(j, nt);
          if (value != -1)
          {
            v += log(value);
          }
          else
          {
            v += log(0.001);
          }
        }
        else
        {
          // this is a non-ACGT symbol
          // => score with low value
          v += log(0.001);
        }
      }
      // cerr << v << " ";
      *(fIter) = v;
      fIter++;
    }
    // cerr << endl;
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
  void SegmentLogBackgroundScore(ScoreVector & logPbx,
                                 SequenceObject * pSeq,
                                 strand_modes s,
                                 BackgroundModel * pBgM, int W) {

    // variables
    int nt,
      content = 0;
    double sum = 0,
      value = 0;

    // model parameters
    int order = pBgM->GetOrder();
    int length = pSeq->Length();

    // check the parameters
    if (order >= length || length < W || W < order)
    {
      // so far we do not do anything with the logPbx
      return;
    }
    int Lw = length - W + 1;

    // iterator
    vector < double >::iterator fIter = logPbx.begin();
    // cerr << "logPbx = ";
    for (int i = 0; i < order; i++)
    {
      nt = pSeq->GetNucleotideAt(s, i);

      // reset sum
      sum = 0;
      for (int j = i; j < order; j++)
      {
        nt = pSeq->GetNucleotideAt(s, j);
        if (nt != -1)
        {
          value = pBgM->GetSnfValueAt(nt);
          if (value != -1)
          {
            sum += log(value);
          }
          else
          {
            sum += log(0.25);
          }
        }
        else
        {
          sum += log(0.25);
        }
      }

      // use mtrans for order -> W + i
      for (int j = order; j < W+i; j++)
      {
        nt = pSeq->GetNucleotideAt(s, j);
        if (nt != -1)
        {
          content = Site2Index(pSeq, s, j - order, order);
          if (content != -1)
          {
            value = pBgM->GetTransitionMatrixValueAt(content, nt);
            if (value != -1)
            {
              sum += log(value);
            }
            else
            {
              value = pBgM->GetSnfValueAt(nt);
              if (value != -1)
              {
                sum += log(value);
              }
              else
              {
                sum += log(0.25);
              }
            }
          }
          else
          {
            sum += log(0.25);
          }
        }
        else
        {
          sum += log(0.25);
        }
      }
      // update value in pbx
      // cerr << sum << " ";
      *(fIter + i) = sum;
    }
    
    // second part
    for (int i = order; i < length; i++)
    {
      nt = pSeq->GetNucleotideAt(s, i);
      if (i < Lw)
      {
        sum = 0;
        for (int j = 0; j < W; j++)
        {
          nt = pSeq->GetNucleotideAt(s, i + j);
          if (nt != -1)
          {
            content = Site2Index(pSeq, s, i + j - order, order);
            if (content != -1)
            {
              value = pBgM->GetTransitionMatrixValueAt(content, nt);
              if (value != -1)
              {
                sum += log(value);
              }
              else
              {
                value = pBgM->GetSnfValueAt(nt);
                sum += log(value);
              }
            }
            else
            {
              sum += log(0.25);
            }
          }
          else
          {
            sum += log(0.25);
          }
        }
        // cerr << sum << " ";
        *(fIter + i) = sum;
      }
      else
      {
        // cerr << endl << "--Warning-- INCLUSIVE::SegmentLogBackgroundScore() outside computation range." << endl;
        *(fIter + i) = W * log(0.25);
        // this to
      }
    }
    // cerr << endl;
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
    SequenceLogBackgroundScore(SequenceObject * pSeq,
                               strand_modes s, BackgroundModel * pBgM) {
    int content,
      nt;
    double logP0 = 0,
      value = 0;
    int order = pBgM->GetOrder();
    int length = pSeq->Length();
    for (int i = 0; i < order; i++)
    {
      nt = pSeq->GetNucleotideAt(s, i);
      if (nt != -1)
      {
        value = pBgM->GetSnfValueAt(nt);
        if (value != -1)
        {
          logP0 += log(value);
        }
        else
        {
          logP0 += log(0.25);
        }
      }
      else
      {
        logP0 += log(0.25);
      }
    }
    for (int i = order; i < length; i++)
    {
      nt = pSeq->GetNucleotideAt(s, i);
      if (nt != -1)
      {
        // get the index of the order previous bases
        content = Site2Index(pSeq, s, i - order, order);
        if (content != -1)
        {
          value = pBgM->GetTransitionMatrixValueAt(content, nt);
          if (value != -1)
          {
            logP0 += log(value);
          }
          else
          {
            logP0 += log(0.25);
          }
        }
        else
        {
          value = pBgM->GetSnfValueAt(nt);
          if (value != -1)
          {
            logP0 += log(value);
          }
          else
          {
            // this is a non-ACGT symbol add 0.25 as background score
            logP0 += log(0.25);
          }
        }
      }
      else
      {
        // this is a non-ACGT symbol add 0.25 as background score
        logP0 += log(0.25);
      }
    }
    return logP0;
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
    ComputeNInstancesProbability(ScoreVector * pExpWx, int n, int L, int W)
  {

    // variables
    int i = 0;
    int ll = L - W + 1;         // array length
    int lw = L - (n * W) + 1;
    if (lw <= 1)
    {
      cerr << "Warning\n CopyProbability: sequence too short: " << L << " < "
        << n << "x" << W << endl;
      return -HUGE_VAL;
    }
    // pointers to store intermediate results
    double sum = 0;
    vector < double >::iterator iterWx = pExpWx->begin();
    vector < double >::iterator iterSum,
      iterPt1,
      iterPt2;

    // number of possible combinations
    double C = log((double) lw);
    for (i = 2; i <= n; i++)
      C += log((double) (L - (n * W) + i)) - log((double) i);

    // define intermediate variables
    vector < double >*pt1 = new vector < double >(ll);
    iterPt1 = pt1->begin();
    vector < double >*pt2 = new vector < double >(ll);
    iterPt2 = pt2->begin();

    // set sum iterator
    iterSum = iterWx + (n - 1) * W;
    for (i = 1; i < n; i++)
    {
      // fill in last element in row
      *(iterPt1 + lw - 1) = *(iterSum + lw - 1);
      *(iterPt2 + lw - 1) = *(iterPt1 + lw - 1) * *(iterWx + (n - i) * W + lw - W);
      for (int j = 1; j <= (lw - 1); j++)
      {
        *(iterPt1 + lw - 1 - j) =
          *(iterPt1 + lw - j) + (*(iterSum + lw - 1 - j));
        *(iterPt2 + lw - 1 - j) =
          *(iterPt1 + lw - 1 - j) * *(iterWx + (n - i) * W + lw - 1 - j - W);
      }
      // let psum point to first element in freshly created row
      iterSum = iterPt1;
    }
    // sum over array of compute scores
    sum = SumArray(iterSum, 0, lw - 1);
    if (sum == 0)
    {
      sum = -HUGE_VAL;
    }
    else
    {
      sum = log(sum) - C;
    }

    // clean up the variables
    delete pt1;
    pt1 = NULL;
    delete pt2;
    pt2 = NULL;
    return sum;
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
    ResetScoreVector(ScoreVector & pVec) 
  {
    for (uint i = 0; i < pVec.size(); i++)
      pVec[i] = 0;
    return;
  }

/******************************************************************************
  Method:       ResetScoreVector
  Class:        
  Arguments:    ScoreVector &pVec
                int value
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
  void
    ResetScoreVector(ScoreVector & pVec, int value) 
  {
    for (uint i = 0; i < pVec.size(); i++)
      pVec[i] = value;
    return;
  }

/******************************************************************************
  Method:       NormalizeScoreVector
  Class:        
  Arguments:    ScoreVector &pVec
  
  Description:  Normalize vector
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
  void
    NormalizeScoreVector(ScoreVector & pVec) 
  {
    double sum = 0;
    // compute sum over array
    for (uint i = 0; i < pVec.size(); i++)
      sum += pVec[i];

    for (uint i = 0; i < pVec.size(); i++)
      pVec[i] /= sum;
    return;
  }

/******************************************************************************
  Method:       ComputePriorDistribution
  Class:        
  Arguments:    ScoreVector &pPrior
                int n
                double prior
  
  Description:  computes a prior distribution of length n starting
                with a given prior. 
                Results are stored in pPrior
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
  void
    ComputePriorDistribution(ScoreVector & pPrior, int n, double prior)
  {
    // reset prior
    ResetScoreVector(pPrior);
    if ((int) pPrior.size() < n + 1)
      pPrior.resize(n + 1, 0);

    // set first two values
    pPrior[0] = 1 - prior;
    pPrior[1] = prior;

    // update next values
    for (int i = 2; i == n; i++)
      pPrior[i] = pPrior[i - 1]/4;

    // normalize prior
    NormalizeScoreVector(pPrior);
    return;
  }


  
}                               // end of the namespace
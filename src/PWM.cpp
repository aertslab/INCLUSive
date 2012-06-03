// 21 july 2009 - 3.1.5

#include "PWM.h"
#include "InstanceMap.h"
#include "Instance.h"

#include <math.h>

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM::PWM(PWM & pCopy)
{
  // motif length
  _length = pCopy.Length();

  // define matrix and initialize with all 0.25
  _pMatrix = new double*[_length];
  for (int i = 0; i < _length; i++)
  {
    _pMatrix[i] = new double[4];
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = 0.25;
  }


  // copy pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = pCopy.GetPseudoCountAt(j);

  // copy matrix
  for (int i = 0; i < _length; i++)
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = pCopy.GetValueAt(i, j);

  _score = 0;
  SetID(pCopy.GetID());
  SetConsensus(pCopy.GetConsensus());

}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM::PWM(int W, Matrix pM)
{
  double sum = 0;
  // motif length
  _length = W;

  // define matrix and initialize to all 0.25
  _pMatrix = new double*[_length];
  for (int i = 0; i < _length; i++)
  {
    _pMatrix[i] = new double[4];
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = 0.25;
  }
  
  for (int i = 0; i < W; i++)
  {
    // calculate sum over row
    sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pM[i][j];
    
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = pM[i][j] / sum;
  }

  // set pseudocounts to very small value
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = 0.0001;

  _score = 0;

  // set identifier and consensus
  _consensus = NULL;
  _ComputeConsensus();
  _id = NULL;

}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM::PWM(int W, Matrix pM, double *snf)
{
  double sum = 0;

  // motif length
  _length = W;

  // define uniform pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  // define matrix
  _pMatrix = new double*[_length];
  for (int i = 0; i < _length; i++)
  {
    _pMatrix[i] = new double[4];
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = 0.25;
  }


  for (int i = 0; i < W; i++)
  {
    // calculate sum over row
    sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pM[i][j] + _pPseudo[j];

    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = (pM[i][j] + _pPseudo[j]) / (sum);
  }

  // set identifier and consensus
  _consensus = NULL;
  _ComputeConsensus();
  _id = NULL;

}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM::PWM(int W, Counts pC)
{
  double sum = 0;
  // motif length
  _length = W;

  // define matrix
  _pMatrix = new double*[_length];
  for (int i = 0; i < _length; i++)
  {
    _pMatrix[i] = new double[4];
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = 0.25;
  }


  // update values
  for (int i = 0; i < W; i++)
  {
    // calculate sum over row
    sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j];

    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = pC[i][j] / sum;
  }

  // set pseudocounts to very small value
    for (int j = 0; j < 4; j++)
    _pPseudo[j] = 0.000001;

  _score = 0;

  // set identifier and consensus
  _consensus = NULL;
  _ComputeConsensus();
  _id = NULL;

}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM::PWM(int W, Counts pC, double *snf)
{
  double sum = 0;
  // motif length
  _length = W;

  // define pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  // define matrix
  _pMatrix = new double*[_length];
  for (int i = 0; i < _length; i++)
  {
    _pMatrix[i] = new double[4];
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = 0.25;
  }

  // update values
  for (int i = 0; i < W; i++)
  {
    // calculate sum over row
    sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j] + _pPseudo[j];

    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = (pC[i][j] + _pPseudo[j]) / (sum);
  }

  _score = 0;

  // set identifier and consensus
  _consensus = NULL;
  _ComputeConsensus();
  _id = NULL;

}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM::PWM(InstanceMap * pMap, double *snf, int w)
{
  double sum = 0;

  _length = w;

  // define pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  // create empty motif count
  Counts pC = new int*[_length];
  for (int i = 0; i < _length; i++)
  {
    pC[i] = new int[4];
    for (int j = 0; j < 4; j++)
      pC[i][j] = 0;
  }


  // get initial element in list
  list < Instance * >::iterator iter = pMap->begin();
  int nt;

  while (iter != pMap->end())
  {
    // get site
    if ((*iter)->Site() != NULL)
    {
      for (int j = 0; j < _length; j++)
      {
        nt = (*((*iter)->Site()))[j];
        if (nt >= 0 && nt < 4)
          pC[j][nt] += 1;
      }
    }
    iter++;                     // next
  }

  // define matrix
  _pMatrix = new double*[_length];
  for (int i = 0; i < _length; i++)
  {
    _pMatrix[i] = new double[4];
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = 0.25;
  }

  // update values
  for (int i = 0; i < _length; i++)
  {
    // calculate sum over row
    sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j] + snf[j];

    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = (pC[i][j] + snf[j]) / sum;

  }

  for (int i = 0; i < _length; i++)
    delete[] pC[i];
  delete[] pC;

  _consensus = NULL;
  _ComputeConsensus();
  _id = NULL;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM::~PWM()
{
  // delete matrix row-by-row
  if (_pMatrix != NULL)
  {
    for (int i=0; i<_length; i++)
      delete[] _pMatrix[i];
    delete[] _pMatrix;
  }
  _pMatrix = NULL;

  // delete identifiers
  if (_id != NULL)
    delete _id;
  _id = NULL;

  if (_consensus != NULL)
    delete _consensus;
  _consensus = NULL;

/*   if ( _pPseudo != NULL )
    delete[] _pPseudo;
 */	
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM *
PWM::RebuildMatrix(int W, Counts pC, double *snf)
{
  double sum = 0;
  if (W != _length)
  {
    // delete old matrix
    for (int i = 0; i < _length; i++)
      delete[] _pMatrix[i];
    delete[]_pMatrix;

    _length = W;
  
    // define new matrix
    _pMatrix = new double*[_length];
    for (int i = 0; i < _length; i++)
    {
      _pMatrix[i] = new double[4];
      for (int j = 0; j < 4; j++)
        _pMatrix[i][j] = 0.25;
    }
  }

  // reset pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  for (int i = 0; i < W; i++)
  {
    // recalculate sum over row
    sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j] + _pPseudo[j];

    // reset matrix elements
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = (pC[i][j] + _pPseudo[j]) / (sum);

  }

  _ComputeConsensus();

  return this;
}

PWM *
PWM::RebuildMatrix(int W, Matrix pC, double *snf)
{
  double sum = 0;
  if (W != _length)
  {
    // delete old matrix
    for (int i = 0; i < _length; i++)
      delete[] _pMatrix[i];
    delete[]_pMatrix;

    _length = W;
  
    // define new matrix
    _pMatrix = new double*[_length];
    for (int i = 0; i < _length; i++)
    {
      _pMatrix[i] = new double[4];
      for (int j = 0; j < 4; j++)
        _pMatrix[i][j] = 0.25;
    }
  }

  // reset pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  for (int i = 0; i < W; i++)
  {
    // recalculate sum over row
    sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j] + _pPseudo[j];

    // reset matrix elements
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = (pC[i][j] + _pPseudo[j]) / (sum);

  }

  _ComputeConsensus();

  return this;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM *
PWM::RebuildMatrix(int W, Counts pC)
{
  double sum = 0;

  if (W != _length)
  {
    // delete old matrix
    for (int i = 0; i < _length; i++)
      delete[] _pMatrix[i];
    delete[]_pMatrix;

    _length = W;
  
    // define new matrix
    _pMatrix = new double*[_length];
    for (int i = 0; i < _length; i++)
    {
      _pMatrix[i] = new double[4];
      for (int j = 0; j < 4; j++)
        _pMatrix[i][j] = 0.25;
    }
  }

  for (int i = 0; i < W; i++)
  {
    // recalculate sum over row
    sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j];

    // reset matrix elements
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = pC[i][j] / sum;

  }

  _ComputeConsensus();

  return this;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM *
PWM::RebuildMatrix(int W, Matrix pC)
{
  if (W != _length)
  {
    // delete old matrix
    for (int i = 0; i < _length; i++)
      delete[] _pMatrix[i];
    delete[]_pMatrix;

    _length = W;
  
    // define new matrix
    _pMatrix = new double*[_length];
    for (int i = 0; i < _length; i++)
    {
      _pMatrix[i] = new double[4];
      for (int j = 0; j < 4; j++)
      {
        _pMatrix[i][j] = 0.25;
      }
    }
  }

  for (int i = 0; i < W; i++)
  {
    // recalculate sum over row
    double sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j];

    // reset matrix elements
    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = pC[i][j] / sum;

  }

  _ComputeConsensus();

  return this;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
string *
PWM::GetConsensus()
{
  if (_consensus == NULL)
    _ComputeConsensus();

  return _consensus;
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
  PWM::GetValueAt(int i, int j)
{
  if (i >= 0 && i < _length && j >= 0 && j < 4)
  {
    return _pMatrix[i][j];
  }
  else
  {
    cerr <<
      "Error\n PWM::GetValueAt() Trying to access index out of range: index = "
      << i << "," << j << endl;
    return -1;
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
  PWM::GetPseudoCountAt(int i)
{
  if (i >= 0 && i < 4)
  {
    return _pPseudo[i];
  }
  else
  {
    cerr <<
      "Error\n PWM::GetPseudoCountAt() Trying to access index out of range: index = "
      << i << endl;
    return -1;
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
PWM *
PWM::SetConsensus(string * cons)
{
  // cerr << "DEBUG PWM::SetConsensus(): --" << *cons << "--" << endl;
  if ( cons->size() == (uint)_length )
  {
    if (_consensus != NULL)
      delete _consensus;
    _consensus = new string(*cons);
  }
  else
  {
    _ComputeConsensus();
  }
  return this;
}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM *
PWM::SetID(string * id)
{
  if (_id != NULL)
    delete _id;

  _id = new string(*id);

  return this;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWM *
PWM::SetScore(double sc)
{
  _score = sc;
  return this;
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
  PWM::_ComputeConsensus()
{
  // reset consensus string
  if (_consensus != NULL)
    delete _consensus;

  _consensus = new string(_length, 'n');

  // variables to store sorted values
  int
    i1 = 0,
    i2 = 0;
  double
    m1 = 0,
    m2 = 0;
  char a[5] = "ACGT";
  char d[8] = "mrwsy.k";

  for (int i = 0; i < _length; i++)
  {
    m1 = 0;
    m2 = 0;
    // find 2 highest scoring bases
    for (int j = 0; j < 4; j++)
    {
      if (_pMatrix[i][j] > m1)
      {
        m2 = m1;                // previous highest scoring becomes second
        i2 = i1;
        m1 = _pMatrix[i][j];    // new highest scoring
        i1 = j;
      }
      else if (_pMatrix[i][j] > m2)
      {
        m2 = _pMatrix[i][j];
        i2 = j;
      }
    }
    int
      index = 0;
    // add letter to string
    if (m1 > 0.66)
    {
      _consensus->replace(i, 1, 1, a[i1]);
    }
    else if (m2 > 0.33)
    {
      if (i1 < i2)
      {
        index = 2 * i1 + i2 - 1;
      }
      else
      {
        index = 2 * i2 + i1 - 1;
      }
      _consensus->replace(i, 1, 1, d[index]);
    }
    else
    {
      _consensus->replace(i, 1, 1, 'n');
    }
  }
  // cerr << "New consensus: " << *_consensus << endl; 
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
  PWM::StderrPrintMatrix()
{
  if (_pMatrix == NULL)
  {
    cerr << "--Error-- PWM::PrintMatrix(): empty matrix" << endl;
    return;
  }
  cerr << "--------------------------------" << endl;
  for (int j = 0; j < 4; j++)
  {
    for (int i = 0; i < _length; i++)
    {
      cerr << _pMatrix[i][j] << "\t";
    }
    cerr << endl;
  }
  cerr << "--------------------------------" << endl;

  return;
}

double
  PWM::ConsensusScore()
{
  double
    score = 0;

  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
      score += _pMatrix[i][j] * log(_pMatrix[i][j]) / log(2.0);
  }

  score = 2 + (score / _length);

  return score;
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
  PWM::InformationContent(double *snf)
{
  double
    score = 0;
  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
      score += _pMatrix[i][j] * (log(_pMatrix[i][j]) -
                                 log(snf[j])) / log(2.0);
  }

  score = score / _length;

  return score;
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
  PWM::MaxScore()
{
  double
    score = 0;
  for (int i = 0; i < _length; i++)
  {
    double
      elem = 0;
    for (int j = 0; j < 4; j++)
    {
      if (_pMatrix[i][j] > elem)
      {
        elem = _pMatrix[i][j];
      }
    }
    score += elem;
  }

  return score;
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
  PWM::MutualInformation(PWM * sbjct, int shift)
{
  int
    leftShift = 0,
    rightShift = 0,
    a1 = 0,
    a2 = 0,
    b1 = 0,
    b2 = 0,
    i = 0,
    j = 0;
  int
    Ws = sbjct->Length();
  double
    mi = 100;
  double
    newMI = 0;

  // 2009/09/03
 _MIshift = 0;

  // compute shift in both directions
  if (_length > Ws)
  {
    leftShift = -shift - _length + Ws;
    rightShift = shift;
  }
  else
  {
    leftShift = -shift;
    rightShift = shift - _length + Ws;
  }

  for (i = leftShift; i <= rightShift; i++)
  {
    // define indices of matching motif parts a1-a2 <-> b1-b2
    a1 = 0 > -i ? 0 : -i;
    a2 = _length < Ws - i ? _length : Ws - i;
    b1 = 0 > i ? 0 : i;
    b2 = Ws < _length + i ? Ws : _length + i;

    newMI = 0;
    for (j = 0; j < a2 - a1; j++)
    {
      newMI += _pMatrix[a1 +
                        j][0] * (log(_pMatrix[a1 + j][0]) -
                                 log(sbjct->
                                     GetValueAt(b1 + j, 0))) / log(2.0);
      newMI += _pMatrix[a1 +
                        j][1] * (log(_pMatrix[a1 + j][1]) -
                                 log(sbjct->
                                     GetValueAt(b1 + j, 1))) / log(2.0);
      newMI += _pMatrix[a1 +
                        j][2] * (log(_pMatrix[a1 + j][2]) -
                                 log(sbjct->
                                     GetValueAt(b1 + j, 2))) / log(2.0);
      newMI += _pMatrix[a1 +
                        j][3] * (log(_pMatrix[a1 + j][3]) -
                                 log(sbjct->
                                     GetValueAt(b1 + j, 3))) / log(2.0);
    }
    // take average over all positions
    newMI = newMI / (a2 - a1);

    // update minimal if new value is smaller
    if (newMI < mi) { mi = newMI; _MIshift = i;}
    //mi = (newMI < mi) ? newMI : mi; 

    // compare with reverse complement
    newMI = 0;
    for (j = 0; j < a2 - a1; j++)
    {
      newMI += _pMatrix[a1 +
                        j][0] * (log(_pMatrix[a1 + j][0]) -
                                 log(sbjct->
                                     GetValueAt(b2 - j - 1, 3))) / log(2.0);
      newMI += _pMatrix[a1 +
                        j][1] * (log(_pMatrix[a1 + j][1]) -
                                 log(sbjct->
                                     GetValueAt(b2 - j - 1, 2))) / log(2.0);
      newMI += _pMatrix[a1 +
                        j][2] * (log(_pMatrix[a1 + j][2]) -
                                 log(sbjct->
                                     GetValueAt(b2 - j - 1, 1))) / log(2.0);
      newMI += _pMatrix[a1 +
                        j][3] * (log(_pMatrix[a1 + j][3]) -
                                 log(sbjct->
                                     GetValueAt(b2 - j - 1, 0))) / log(2.0);
    }
    // take average over all positions
    newMI = newMI / (a2 - a1);

    // update minimal if new value is smaller
    if (newMI < mi) { mi = newMI; _MIshift = i;}
    //mi = (newMI < mi) ? newMI : mi;
  }

  return mi;

}

/******************************************************************************
  Description:  compute BLiC score between two matrices, 
                also return alignshift/overlap/strand/pvalue/common
  Author_Date:  MC_2009/07/21
  Author_Date:  MC_2011/03/28 : update BLiC formula (rewrite to clear format)
******************************************************************************/
double *
PWM::BLiC(PWM *dbMatrix, int shift, double *dirichlet, double *bgsnf, bool a2a1)
{
  int leftShift(0),rightShift(0),a1(0),a2(0),b1(0),b2(0),i(0),j(0);
  int Ws = dbMatrix->Length();
  double BLiC(-100), newBLiC(0);
  double BLiC_shift = 100; //
  double BLiC_overlap = 100; //
  double BLiC_strand = 100; // 
  double dirsum = 0;
  for(int i = 0; i < 4; i++) dirsum += dirichlet[i];
  double weight2 = dbMatrix->GetWeight();
  Matrix dbmatrix = dbMatrix->GetMatrix();
  double countN1, countN2, countN3, freq3, freq1, freq2;// '3' = common PWM
  double weight3 = _weight + weight2;// '3' = common PWM
	
	
  // compute shift in both directions //klopt
  if (_length > Ws) // L1 > L2
  {
    leftShift = -shift - _length + Ws;
    rightShift = shift;
  }
  else
  {
    leftShift = -shift;
    rightShift = shift - _length + Ws;
  }

  for (i = leftShift; i <= rightShift; i++)
  {
    // define indices of matching motif parts a1-a2 <-> b1-b2
    a1 = 0 > -i ? 0 : -i; // if (i > 0) a1 = 0; else a1 = -i;
    a2 = _length < Ws - i ? _length : Ws - i;
          //if ( L2 -i > L1) a2 = L1; else a2 = L2-i;
    b1 = 0 > i ? 0 : i; //if (i < 0) b1 = 0; else b1 = i; 
    b2 = Ws < _length + i ? Ws : _length + i;
          // if (L1 + i > L2) b2 = L2; else b2 = L1 + i; 

    newBLiC = 0;
    for (j = 0; j < a2 - a1; j++)
    {     
      for(int n = 0; n < 4; n++) // n = a,c,g,t
      {
        countN1 = _weight * _pMatrix[a1+j][n];
        countN2 = weight2 * dbmatrix[b1+j][n];
        //newBLiC += 2*(countN1+countN2)*(log(countN1+countN2+dirichlet[n])-
         //                              log(_weight + weight2 + dirsum))-
          //      countN1*(log(countN1+dirichlet[n])-log(_weight+dirsum))-
        //        countN2*(log(countN2+dirichlet[n])-log(weight2+dirsum))-
        //        (countN1+countN2)*log(bgsnf[n]);
		  
// format change : rewrite in more clear format : see Guidelines fig 6

        countN3 = countN1 + countN2;
        freq3 = (countN3 + dirichlet[n])/(weight3 + dirsum);
		freq1 = (countN1 + dirichlet[n])/(_weight + dirsum);
		freq2 = (countN2 + dirichlet[n])/(weight2 + dirsum);

        newBLiC += countN1*log(freq3/freq1) 
			       + countN2*log(freq3/freq2)
			       + countN3*log(freq3/bgsnf[n]);
      }
    }
    // take average over all positions
    if (a2a1) {newBLiC = newBLiC / (a2 - a1);}

    // update maximum if new value is higher
    if (newBLiC > BLiC) 
    { BLiC = newBLiC; BLiC_shift = i; BLiC_overlap = (a2-a1); BLiC_strand = 1;} 

    // compare with reverse complement
    newBLiC = 0;
    for (j = 0; j < a2 - a1; j++)
    {     
     for(int n = 0; n < 4; n++) // n = a,c,g,t
     {
       countN1 = _weight * _pMatrix[a1+j][n];
       countN2 = weight2 * dbmatrix[b2-j-1][3-n];
      
       //newBLiC += 2*(countN1+countN2)*(log(countN1+countN2+dirichlet[n])-
        //                               log(_weight + weight2 + dirsum))-
        //        countN1*(log(countN1+dirichlet[n])-log(_weight+dirsum))-
       //         countN2*(log(countN2+dirichlet[n])-log(weight2+dirsum))-
       //         (countN1+countN2)*log(bgsnf[n]);
// change exactly same as above 
        countN3 = countN1 + countN2;
        freq3 = (countN3 + dirichlet[n])/(weight3 + dirsum);
		freq1 = (countN1 + dirichlet[n])/(_weight + dirsum);
		freq2 = (countN2 + dirichlet[n])/(weight2 + dirsum);

        newBLiC += countN1*log(freq3/freq1) 
			       + countN2*log(freq3/freq2)
			       + countN3*log(freq3/bgsnf[n]);
      }
    }
    // take average over all positions (default is set false)
    if (a2a1) {newBLiC = newBLiC / (a2 - a1);}

    // update maximum if new value is higher
    if (newBLiC > BLiC) 
    { BLiC = newBLiC; BLiC_shift = i; BLiC_overlap = (a2-a1); BLiC_strand = -1;} 
  }

  double * output = new double[4]; 
  output[0] = BLiC; output[1] = BLiC_shift; output[2] = BLiC_overlap; 
  output[3]= BLiC_strand;
  return output;
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
  PWM::WeightedKullbackLeiberDistance(PWM* sbjct, vector< double >* pVec)
{
  double
    newMI = 0,
    mi = 0,
    sum = 0;
  
  for (int j = 0; j < _length; j++)
  {
    newMI = 0;
    newMI += _pMatrix[j][0] * (log(_pMatrix[j][0]) -
                               log(sbjct->GetValueAt(j, 0))) / log(2.0);
    newMI += _pMatrix[j][1] * (log(_pMatrix[j][1]) -
                               log(sbjct->GetValueAt(j, 1))) / log(2.0);
    newMI += _pMatrix[j][2] * (log(_pMatrix[j][2]) -
                               log(sbjct->GetValueAt(j, 2))) / log(2.0);
    newMI += _pMatrix[j][3] * (log(_pMatrix[j][3]) -
                               log(sbjct->GetValueAt(j, 3))) / log(2.0);

    mi += newMI * (*pVec)[j];
    
    sum += (*pVec)[j];
    
  }

  // normalize mi with sum of weights
  mi = mi / sum;
  
  return mi;

}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/

char 
PWM::GetConsensusSymbolAt(int index)
{
  if ( index < 0 || index > _length )
  {
    return 'n';
  }
  if ( _consensus == NULL || _consensus->size() != (uint)_length)
    _ComputeConsensus();
  
  return _consensus->at(index);
}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/

PWM*
PWM::SubMatrix(int index, int length)
{
  PWM * pwm = NULL;
  
  if ( index < 0 || index > _length || length < 0 || length > _length || index + length > _length )
  {
    cerr << "Warning PWM::SubMatrix: Unable to create submatrix, index out of bounds." << endl;
    return NULL;
  }    
    
  // define sub matrix
  Matrix temp = new double*[length];
  for (int i = 0; i < length; i++)
  {
    temp[i] = new double[4];
    for (int j = 0; j < 4; j++)
      temp[i][j] = _pMatrix[index+i][j];
  }
  
  // create new matrix
  pwm = new PWM(length, temp);

  for (int i = 0; i < length; i++)
    delete[] temp[i];
  delete[] temp;

  return pwm;
}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/


PWM*
PWM::ReverseSubMatrix(int index, int length)
{
  PWM * pwm = NULL;
  
  if ( index > _length || length < 0 || length > _length || index - length + 1 < 0 )
  {
    cerr << "Warning PWM::ReverseSubMatrix: Unable to create submatrix, index out of bounds." << endl;
    return NULL;
  }    
    
  // define sub matrix
  Matrix temp = new double*[length];
  for (int i = 0; i < length; i++)
  {
    temp[i] = new double[4];
    for (int j = 0; j < 4; j++)
      temp[i][j] = _pMatrix[index-i][3-j];
  }
  
  // create new matrix object
  pwm = new PWM(length, temp);

  for (int i = 0; i < length; i++)
    delete[] temp[i];
  delete[] temp;
  
  return pwm;
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
PWM::ShuffleMatrix()
{
  int value = 0;
  double row[4];
  
  // shuffle matrix rows
  for (int i = 0; i < _length; i++)
  {
    // switch row i with row 'value'
    value = (int)(_length * (random()/(RAND_MAX + 1.0)));
    for (int j = 0; j < 4; j++)
    {
      row[j] = _pMatrix[i][j];   // temporary copy of row i
      _pMatrix[i][j] = _pMatrix[value][j];
      _pMatrix[value][j] = row[j]; 
    }
  }
  // add 03/09/2009
  _ComputeConsensus();
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


PWM*
PWM::NewShuffledMatrix()
{

  int indices[_length];
  int value = 0, dummy = 0;
  
  // create initial indices
  for ( int i=0; i<_length; i++)
    indices[i] = i;

  // shuffle indices
  for ( int i=0; i<_length; i++ )
  {
    value = (int)(_length * (random()/(RAND_MAX + 1.0)));
    dummy = indices[value];
    indices[value] = indices[i];
    indices[i] = dummy;
  }

  // shuffle matrix rows
  Matrix temp = new double*[_length];
  for (int i = 0; i < _length; i++)
  {
    temp[i] = new double[4];
    for (int j = 0; j < 4; j++)
      temp[i][j] = _pMatrix[indices[i]][j];
  }
  
  // create new motif model from shuffled matrix
  PWM * pwm = new PWM(_length, temp);

  // cleanup local variable
  for (int i = 0; i < _length; i++)
    delete[] temp[i];
  delete[] temp;


  return pwm;
}

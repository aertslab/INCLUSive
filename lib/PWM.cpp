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
  _pMatrix = new double[_length][4];
  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = 0.25;
    }
  }


  // copy pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = pCopy.GetPseudoCountAt(j);

  // copy matrix
  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = pCopy.GetValueAt(i, j);
    }
  }

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
  // motif length
  _length = W;

  // define matrix and initialize to all 0.25
  _pMatrix = new double[W][4];
  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = 0.25;
    }
  }

  for (int i = 0; i < W; i++)
  {
    // calculate sum over row
    double sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pM[i][j];

    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = pM[i][j] / sum;
    }
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
PWM::PWM(int W, Matrix pM, double *snf)
{
  // motif length
  _length = W;

  // define uniform pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  // define matrix
  _pMatrix = new double[W][4];
  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = 0.25;
    }
  }


  for (int i = 0; i < W; i++)
  {
    // calculate sum over row
    double sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pM[i][j] + _pPseudo[j];

    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = (pM[i][j] + _pPseudo[j]) / (sum);
    }
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
  // motif length
  _length = W;

  // define matrix
  _pMatrix = new double[W][4];
  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = 0.25;
    }
  }


  // update values
  for (int i = 0; i < W; i++)
  {
    // calculate sum over row
    double sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j];

    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = pC[i][j] / sum;
    }
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
  // motif length
  _length = W;

  // define pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  // define matrix
  _pMatrix = new double[W][4];
  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = 0.25;
    }
  }

  // update values
  for (int i = 0; i < W; i++)
  {
    // calculate sum over row
    double sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j] + _pPseudo[j];

    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = (pC[i][j] + _pPseudo[j]) / (sum);
    }
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
  // set identifier and consensus
  _consensus = NULL;
  _id = NULL;
  _length = w;

  // define pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  // create empty motif count
  Counts pC = new int[_length][4];
  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      pC[i][j] = 0;
    }
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
  _pMatrix = new double[_length][4];
  for (int i = 0; i < _length; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      _pMatrix[i][j] = 0.25;
    }
  }

  // update values
  for (int i = 0; i < _length; i++)
  {
    // calculate sum over row
    double sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pC[i][j] + snf[j];

    for (int j = 0; j < 4; j++)
      _pMatrix[i][j] = (pC[i][j] + snf[j]) / sum;

  }

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
  // delete matrices
  // delete[] _pPseudo;
  if (_pMatrix != NULL)
    delete[]_pMatrix;
  _pMatrix = NULL;

  // delete identifiers
  if (_id != NULL)
    delete _id;
  _id = NULL;

  if (_consensus != NULL)
    delete _consensus;
  _consensus = NULL;

  if ( _pPseudo != NULL )
    delete[] _pPseudo;
  
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
  if (W != _length)
  {
    // need to define new matrix
    delete[]_pMatrix;
    _pMatrix = new double[W][4];

  }

  _length = W;

  // reset pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  for (int i = 0; i < W; i++)
  {
    // recalculate sum over row
    double sum = 0;
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
  if (W != _length)
  {
    // need to define new matrix
    delete[]_pMatrix;
    _pMatrix = new double[W][4];
  }

  _length = W;

  // reset pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = snf[j];

  for (int i = 0; i < W; i++)
  {
    // recalculate sum over row
    double sum = 0;
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
  if (W != _length)
  {
    // need to define new matrix
    delete[]_pMatrix;
    _pMatrix = new double[W][4];

  }

  _length = W;

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
PWM *
PWM::RebuildMatrix(int W, Matrix pC)
{
  if (W != _length)
  {
    // need to define new matrix
    delete[]_pMatrix;
    _pMatrix = new double[W][4];

  }

  _length = W;

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
  char *
    a = "ACGT";
  char *
    d = "mrwsy.k";

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
    mi = (newMI < mi) ? newMI : mi;

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
                                 log(sbjct->GetValueAt(b2 - j - 1, 0))) /
        log(2.0);
    }
    // take average over all positions
    newMI = newMI / (a2 - a1);

    // update minimal if new value is smaller
    mi = (newMI < mi) ? newMI : mi;
  }

  return mi;

}


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


PWM*
PWM::SubMatrix(int index, int length)
{
  PWM * pwm = NULL;
  
  if ( index < 0 || length < 0 || length > _length || index + length > _length )
  {
    cerr << "Warning PWM::SubMatrix: Unable to create submatrix, index out of bounds." << endl;
    return NULL;
  }    
    
  // define sub matrix
  Matrix pMatrix = new double[length][4];
  for (int i = 0; i < length; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      pMatrix[i][j] = _pMatrix[index+i][j];
    }
  }
  
  // create new matrix
  pwm = new PWM(length, pMatrix);
  // cerr << "DEBUG: PWM::SubMatrix(): New matrix: " << *(pwm->GetConsensus()) << endl;

  return pwm;
}

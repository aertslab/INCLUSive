#include "SpacerPWM.h"
#include <math.h>


SpacerPWM::SpacerPWM(int W1, Matrix pM1, int spacer, int W2, Matrix pM2)
{

  // define motif lengths
  _length1 = W1;
  _length2 = W2;

  // define spacer
  _spacer = spacer;

  // define pseudo counts
  for (int j = 0; j < 4; j++)
    _pPseudo[j] = 0.05;

  // define matrix part 1
  _pMatrix1 = new double[W1][4];
  for (int i = 0; i < W1; i++)
  {
    // calculate sum over row
    double sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pM1[i][j] + _pPseudo[j];

    for (int j = 0; j < 4; j++)
      _pMatrix1[i][j] = (pM1[i][j] + _pPseudo[j]) / (sum);
  }

  // define matrix part 2
  _pMatrix2 = new double[W2][4];
  for (int i = 0; i < W1; i++)
  {
    // calculate sum over row
    double sum = 0;
    for (int j = 0; j < 4; j++)
      sum += pM2[i][j] + _pPseudo[j];

    for (int j = 0; j < 4; j++)
      _pMatrix2[i][j] = (pM2[i][j] + _pPseudo[j]) / (sum);
  }

  // set some initial values
  _score = 0;
  _consensus = NULL;
  _id = NULL;


  return;
}


SpacerPWM::~SpacerPWM()
{

  // delete matrix part 1
  if (_pMatrix1 != NULL)
    delete[]_pMatrix1;
  _pMatrix1 = NULL;

  // delete matrix part 2
  if (_pMatrix2 != NULL)
    delete[]_pMatrix2;
  _pMatrix2 = NULL;

  // delete defined strings

}

SpacerPWM *
SpacerPWM::SetConsensus(string * cons)
{
  if (_consensus != NULL)
    delete _consensus;

  _consensus = new string(*cons);
  return this;
}

SpacerPWM *
SpacerPWM::SetID(string * id)
{
  if (_id != NULL)
    delete _id;

  _id = new string(*id);

  return this;
}


SpacerPWM *
SpacerPWM::SetScore(double sc)
{
  _score = sc;
  return this;
}


// private functions to compute consensus
// consensus should have the form M1-nX-M2
void
SpacerPWM::_ComputeConsensus()
{
  //if ( _consensus != NULL )
  //  delete _consensus;
  int totalLength = _length1 + _spacer + _length2;
  char *a = "ACGT";
  char *d = "mrwsy.k";

  //_consensus = new string(_length,'n');
  if (_consensus == NULL)
  {
    _consensus = new string(totalLength, 'n');
  }
  else if (totalLength != (int) _consensus->size())
  {
    _consensus->resize(totalLength);
  }

  // variables to store sorted values
  int i1 = 0, i2 = 0;
  double m1 = 0, m2 = 0;

  for (int i = 0; i < _length1; i++)
  {
    m1 = 0;
    m2 = 0;
    // find 2 highest scoring bases
    for (int j = 0; j < 4; j++)
    {
      if (_pMatrix1[i][j] > m1)
      {
        m2 = m1;                // previous highest scoring becomes second
        i2 = i1;
        m1 = _pMatrix1[i][j];   // new highest scoring
        i1 = j;
      }
      else if (_pMatrix1[i][j] > m2)
      {
        m2 = _pMatrix1[i][j];
        i2 = j;
      }
    }
    int index = 0;
    // add letter to string
    if (m1 > 0.65)
    {
      _consensus->replace(i, 1, a[i1]);
    }
    else if (m2 > 0.35)
    {
      if (i1 < i2)
      {
        index = 2 * i1 + i2 - 1;
      }
      else
      {
        index = 2 * i2 + i1 - 1;
      }
      _consensus->replace(i, 1, d[index]);
    }
    else
    {
      _consensus->replace(i, 1, 'n');
    }
  }

  // part two
  for (int i = 0; i < _length2; i++)
  {
    m1 = 0;
    m2 = 0;
    // find 2 highest scoring bases
    for (int j = 0; j < 4; j++)
    {
      if (_pMatrix2[i][j] > m1)
      {
        m2 = m1;                // previous highest scoring becomes second
        i2 = i1;
        m1 = _pMatrix2[i][j];   // new highest scoring
        i1 = j;
      }
      else if (_pMatrix2[i][j] > m2)
      {
        m2 = _pMatrix2[i][j];
        i2 = j;
      }
    }
    int index = 0;
    // add letter to string
    if (m1 > 0.65)
    {
      _consensus->replace(_length1 + _spacer + i, 1, a[i1]);
    }
    else if (m2 > 0.35)
    {
      if (i1 < i2)
      {
        index = 2 * i1 + i2 - 1;
      }
      else
      {
        index = 2 * i2 + i1 - 1;
      }
      _consensus->replace(_length1 + _spacer + i, 1, d[index]);
    }
    else
    {
      _consensus->replace(_length1 + _spacer + i, 1, 'n');
    }
  }

  // cerr << "New consensus: " << *_consensus << endl;
  return;
}

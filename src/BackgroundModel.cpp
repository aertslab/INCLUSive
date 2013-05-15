#include "BackgroundModel.h"

#include <iostream>
#include <math.h>


// Constructors

/******************************************************************************
  Method:       BackgroundModel
  Class:        BackgroundModel
  Arguments:    int order
                double (* pTrans)[4]
                double *pFreq
                double *pSnf
  Description:  class constructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
BackgroundModel::BackgroundModel(int order, double **pTrans,
                                 double *pFreq, double *pSnf)
{
  int i, j;

  if (order > 0)
  {
    _length = (int) pow(4.0, order);
  }
  else if (order == 0)
  {
    _length = 4;
  }
  else
  { // define this error during input of CreatBackgroundModel. 
    cerr <<
      "--Error-- BackgroundModel::BackgroundModel -> order should be greater than 0"
      << endl;
    order = -1;
    _length = 0;
    _snf = NULL;
    _transitionMatrix = NULL;
    _oligoFrequencyMatrix = NULL;
    return;
  }

  // set order
  _order = order;

  // make copy of snf
  _snf = new double[4];
  for (i = 0; i < 4; i++)
    _snf[i] = pSnf[i];

  // make copy of oligo frequencies
  _oligoFrequencyMatrix = new double[_length];
  for (i = 0; i < _length; i++)
    _oligoFrequencyMatrix[i] = pFreq[i];


  // copy transition matrix
  _transitionMatrix = new double*[_length];
  for (i = 0; i < _length; i++)
  {
    _transitionMatrix[i] = new double[4];
    for (j = 0; j < 4; j++)
      _transitionMatrix[i][j] = pTrans[i][j];
  }

  // set default values for organism
  _organism = NULL;
  _pFile = NULL;
}


/******************************************************************************
  Method:       ~BackgroundModel
  Class:        BackgroundModel
  Arguments:    none
  
  Description:  class destructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
BackgroundModel::~BackgroundModel()
{
  if (_transitionMatrix != NULL)
  {
    for (int i = 0; i < _length; i++)
      delete[] _transitionMatrix[i];
    delete[] _transitionMatrix;
  }
  _transitionMatrix = NULL;

  if (_oligoFrequencyMatrix != NULL)
    delete[] _oligoFrequencyMatrix;
  _oligoFrequencyMatrix = NULL;

  if (_snf != NULL)
    delete[] _snf;
  _snf = NULL;

  
  if ( _pFile != NULL )
    delete _pFile;
  _pFile = NULL;
  if ( _organism != NULL )
    delete _organism;
  _organism = NULL;

}


/******************************************************************************
  Method:       GetSnfValueAt
  Class:        BackgroundModel
  Arguments:    int i
  
  Description:  returns single nucleotide frequency
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
BackgroundModel::GetSnfValueAt(int i)
{
  if (i >= 0 && i < 4)
  {
    return _snf[i];
  }
  else
  {
    return -1;
  }
}

/******************************************************************************
  Method:       GetTransitionMatrixValueAt
  Class:        BackgroundModel
  Arguments:    int i
                int j 
  
  Description:  returns value at position (i,j) from the transition matrix
                of the background model
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
BackgroundModel::GetTransitionMatrixValueAt(int i, int j)
{
  if (i >= 0 && i < _length && j >= 0 && j < 4)
  {
    return _transitionMatrix[i][j];
  }
  else
  {
    return -1;
  }
}
/******************************************************************************
  Description:  return max transition frequency for given nucleotide 
  Author_Date:  MC_2013/02/04
******************************************************************************/
double 
BackgroundModel::GetMaxTransitionMatrixValueAt(int j)
{
  double max = 0.0000000000000000000000000000001;
  for (int i = 0; i < _length; i++)
  {
    if (_transitionMatrix[i][j] > max)
      max = _transitionMatrix[i][j];
  }
  return max;
}
/******************************************************************************
  Description:  return min transition frequency for given nucleotide 
  Author_Date:  MC_2013/02/04
******************************************************************************/
double 
BackgroundModel::GetMinTransitionMatrixValueAt(int j)
{
  double min = 1;
  for (int i = 0; i < _length; i++)
  {
    if (_transitionMatrix[i][j] < min)
      min = _transitionMatrix[i][j];
  }
  return min;
}

/******************************************************************************
  Method:       GetOligoFrequencyValueAt
  Class:        BackgroundModel
  Arguments:    int i
  
  Description:  returns the value at position i in the oligo frequency
                vector
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double
BackgroundModel::GetOligoFrequencyValueAt(int i)
{
  if (i >= 0 && i < _length)
  {
    return _oligoFrequencyMatrix[i];
  }
  else
  {
    return -1;
  }
}

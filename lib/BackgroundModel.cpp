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
BackgroundModel::BackgroundModel(int order, double (*pTrans)[4],
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
  {
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
  _transitionMatrix = new double[_length][4];
  for (i = 0; i < _length; i++)
  {
    for (j = 0; j < 4; j++)
    {
      _transitionMatrix[i][j] = pTrans[i][j];
    }
  }

  // set default values for organism
  _organism = "";
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
    delete[]_transitionMatrix;
  _transitionMatrix = NULL;

  if (_oligoFrequencyMatrix != NULL)
    delete[]_oligoFrequencyMatrix;
  _oligoFrequencyMatrix = NULL;

  if (_snf != NULL)
    delete[]_snf;
  _snf = NULL;

  _pFile = NULL;
}


/******************************************************************************
  Method:       GetOrder
  Class:        BackgroundModel
  Arguments:    none
  
  Description:  returns order of the background model
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
BackgroundModel::GetOrder()
{
  return _order;
}


/******************************************************************************
  Method:       GetOrganism
  Class:        BackgroundModel
  Arguments:    none
  
  Description:  returns organism name
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
string *
BackgroundModel::GetOrganism()
{
  return &_organism;
}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double *
BackgroundModel::GetSNF()
{
  return _snf;
}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
double *
BackgroundModel::GetOligoFrequencyMatrix()
{
  return _oligoFrequencyMatrix;
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

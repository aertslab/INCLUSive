#include "RandomNumber.h"
#include "Distribution.h"
#include <math.h>

/******************************************************************************
  Method: 
  Class:        Distribution
  Arguments:    nt length, float *pArray
  
  Description:  constructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(int length, float *pArray)
{
  if (length > 0)
  {
    // copy distribution
    _dist = new vector<double>(length);

    vector < double >::iterator iter = _dist->begin();
    for (int i = 0; i < length; i++)
    {
      *iter = (double)pArray[i];
      iter++;
    }

    _normalized = _Normalize();

    if (!_normalized)
    {
      // delete vector if not normalized
      cerr << "--Error-- Distribution(): Problems normalizing distribution -> reset to NULL." << endl;
      delete _dist;
      _dist = NULL;
    }

    return;

  }
  else
  {
    _dist = NULL;
    cerr <<
      "--Error-- Distribution(): length of vector should be greater than 0."
      << endl;
  }
}


/******************************************************************************
  Method: 
  Class:        Distribution
  Arguments:    int length, int *pArray
  
  Description:  constructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(int length, int *pArray)
{
  if (length > 0)
  {
    // copy distribution
    _dist = new vector<double>(length);

    vector < double >::iterator iter = _dist->begin();
    for (int i = 0; i < length; i++)
    {
      *iter = (double)pArray[i];
      iter++;
    }

    _normalized = _Normalize();

    if (!_normalized)
    {
      // delete vector if not normalized
      cerr << "--Error-- Distribution(): Problems normalizing distribution -> reset to NULL." << endl;
      delete _dist;
      _dist = NULL;
    }

    return;

  }
  else
  {
    _dist = NULL;
    cerr <<
      "--Error-- Distribution(): length of vector should be greater than 0."
      << endl;
  }
}


/******************************************************************************
  Method: 
  Class:        Distribution
  Arguments:    int length, double *pArray
  
  Description:  constructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(int length, double *pArray)
{
  if (length > 0)
  {
    // copy distribution
    _dist = new vector < double >(length);

    vector < double >::iterator iter = _dist->begin();
    for (int i = 0; i < length; i++)
    {
      *iter = (double) pArray[i];
      iter++;
    }
    _normalized = _Normalize();

    if (!_normalized)
    {
      // delete vector if not normalized
      cerr <<
        "--Error-- Distribution(): Problems normalizing distribution -> reset to NULL."
        << endl;
      delete _dist;
      _dist = NULL;
    }

    return;

  }
  else
  {
    _dist = NULL;
    cerr <<
      "--Error-- Distribution(): length of vector should be greater than 0."
      << endl;
  }
}


/******************************************************************************
  Method: 
  Class:        Distribution
  Arguments:    vector<float> * pArray
  
  Description:  constructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(vector<float> * pArray)
{
  if (pArray->size() > 0)
  {
    // copy distribution
    _dist = new vector < double >;

    // cerr << "vector<float> input: ";
    for (int i = 0; i < (int) pArray->size(); i++)
    {
      _dist->push_back((double) (*pArray)[i]);
      // cerr << (*pArray)[i] << "|" << (*_dist)[i] << " ";
    }
    // cerr << endl;

    // normalize
    _normalized = _Normalize();
    if (!_normalized)
    {
      // delete vector if not normalized
      cerr <<
        "--Error-- Distribution(): Problems normalizing distribution -> reset to NULL."
        << endl;
      delete _dist;
      _dist = NULL;
    }

    return;

  }
  else
  {
    _dist = NULL;
    cerr <<
      "--Error-- Distribution(): length of vector should be greater than 0."
      << endl;
  }
}


/******************************************************************************
  Method: 
  Class:        Distribution
  Arguments:    vector<double> * pArray
  
  Description:  constructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(vector<double> * pArray)
{
  if (pArray->size() > 0)
  {
    // copy distribution
    _dist = new vector < double >;

    // cerr << "vector<double> input: ";
    for (int i = 0; i < (int) pArray->size(); i++)
    {
      _dist->push_back((double) (*pArray)[i]);
      // cerr << (*pArray)[i] << "|" << (*_dist)[i] << " ";
    }
    // cerr << endl;

    _normalized = _Normalize();
    if (!_normalized)
    {
      // delete vector if not normalized
      cerr <<
        "--Error-- Distribution(): Problems normalizing distribution -> reset to NULL."
        << endl;
      delete _dist;
      _dist = NULL;
    }

    return;

  }
  else
  {
    _dist = NULL;
    cerr <<
      "--Error-- Distribution(): length of vector should be greater than 0."
      << endl;
  }
}


/******************************************************************************
  Method: 
  Class:        Distribution         
  Arguments:    vector<int> * pArray
  
  Description:  constructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(vector<int> * pArray)
{
  if (pArray->size() > 0)
  {
    // copy distribution
    _dist = new vector < double >;

    // cerr << "vector<int> input: ";
    for (int i = 0; i < (int) pArray->size(); i++)
    {
      _dist->push_back((double) (*pArray)[i]);
      // cerr << (*pArray)[i] << "|" << (*_dist)[i] << " ";
    }
    // cerr << endl;

    // normalize
    _normalized = _Normalize();
    if (!_normalized)
    {
      // delete vector if not normalized
      cerr <<
        "--Error-- Distribution(): Problems normalizing distribution -> reset to NULL."
        << endl;
      delete _dist;
      _dist = NULL;
    }

    return;

  }
  else
  {
    _dist = NULL;
    cerr <<
      "--Error-- Distribution(): length of vector should be greater than 0."
      << endl;
  }
}


/******************************************************************************
  Method:       
  Class:        Distribution
  Arguments:    none
  
  Description:	destructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::~Distribution()
{
  if (_dist != NULL)
    delete _dist;
  _dist = NULL;
}


/******************************************************************************
  Method:       TakeSample()
  Class:        Distribution
  Arguments:    none
  
  Description:  Take one sample from the given distribution
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
Distribution::TakeSample()
{
  // cerr << "Starting TakeSample()" << endl;
  int sample = -1;

  if (_dist == NULL)
  {
    cerr << "--Warning-- Distribution::TakeSample(): Trying to sample from empty distribution -> return -1."
      << endl;
    return sample;
  }

  // get random number
  double cumsum = 0;
  double value = _rn.GetUniform();
  // cerr << "Sampled value = " << value << endl;

  for (sample = 0; sample < (int) _dist->size(); sample++)
  {
    cumsum += (*_dist)[sample];
    // cerr << sample << ":" << cumsum << "|" << value << " "; 
    if (cumsum >= value)
      break;
  }
  // cerr << endl;

  return sample;
}


/******************************************************************************
  Method:       SelectMax()
  Class:        Distribution
  Arguments:    none
  
  Description:  Get the index of the highest scoring entry
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
Distribution::SelectMax()
{
  int m = 0;
  // return -1 if distribution is empty
  if (_dist == NULL)
    return -1;

  // iterate through list and output index of highest score
  for (int i = 0; i < (int) _dist->size(); i++)
    m = ((*_dist)[m] < (*_dist)[i]) ? i : m;

  return m;
}


/******************************************************************************
  Method:       ApplyMask()
  Class:        Distribution
  Arguments:    int start, int length
  
  Description:  Set the distribution to 0 from position start to start+length-1
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution *
Distribution::ApplyMask(int start, int length)
{
  if (_dist == NULL)
    return NULL;

  // check start position
  if (start < 0)
    start = 0;

  if (start >= (int) _dist->size())     // no need to update mask here
    return this;

  // check end position
  if (start + length > (int) _dist->size())
    length = (int) _dist->size() - start;

  vector < double >::iterator iter = _dist->begin();
  for (int i = 0; i < length; i++)
    *(iter + start + i) = 0;

  // renormalize distribution
  _normalized = _Normalize();
  if (!_normalized)
  {
    // delete vector if not normalized
    cerr << "--Warning-- Distribution::ApplyMask(): all position set to 0 "
         << "-> reset distribution to NULL." << endl;
    // delete _dist;
    // _dist = NULL;
  }

  return this;

}


/******************************************************************************
  Method:       _Normalize()
  Class:        Distribution
  Arguments:    none
  
  Description:  internal function to normalize the distribution
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
bool
Distribution::_Normalize()
{

  // do nothing if distribution is empty
  if (_dist == NULL)
  {
    cerr << "--Warning-- Distribution::_Normalize(): _dist is NULL." << endl;
    return false;
  }

  // compute sum over array
  double sum = 0;
  vector < double >::iterator iter;
  for (iter = _dist->begin(); iter != _dist->end(); iter++)
    sum += *iter;

  if (sum == 0)
    return false;

  // update values
  for (iter = _dist->begin(); iter != _dist->end(); iter++)
    (*iter) = (*iter) / sum;

  return true;
}

/******************************************************************************
  Description:  Distribution(ScoreVector *pArray, int start, int length);
  Author_Date:  MC_2009/09/10 (copy from 5.0.4)
******************************************************************************/
Distribution::Distribution(ScoreVector *pArray, int start, int length)
{
  if ((int)pArray->size() > start)
  {
    // copy distribution
    _dist = new vector < double >;
    for (int i = start; i < start + length; i++) //
      _dist->push_back((*pArray)[i]);
    // normalize
    if (!_Normalize())
    { delete _dist;_dist = NULL; // all zeros
      return;
    }
  }
  else // vector too short
    return;
}

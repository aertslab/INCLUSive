#include "RandomNumber.h"
#include "Distribution.h"
#include <math.h>

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(int length, float *pArray)
{
  if (length > 0)
  {
    // copy distribution
    _dist = new vector < float >(length);

    vector < float >::iterator iter = _dist->begin();
    for (int i = 0; i < length; i++)
    {
      *iter = pArray[i];
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
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(int length, double *pArray)
{
  if (length > 0)
  {
    // copy distribution
    _dist = new vector < float >(length);

    vector < float >::iterator iter = _dist->begin();
    for (int i = 0; i < length; i++)
    {
      *iter = (float) pArray[i];
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
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(vector < float >*pArray)
{
  if (pArray->size() > 0)
  {
    // copy distribution
    _dist = new vector < float >;

    // cerr << "vector<int> input: ";
    for (int i = 0; i < (int) pArray->size(); i++)
    {
      _dist->push_back((float) (*pArray)[i]);
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
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(vector < double >*pArray)
{
  if (pArray->size() > 0)
  {
    // copy distribution
    _dist = new vector < float >;

    // cerr << "vector<int> input: ";
    for (int i = 0; i < (int) pArray->size(); i++)
    {
      _dist->push_back((float) (*pArray)[i]);
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
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Distribution::Distribution(vector < int >*pArray)
{
  if (pArray->size() > 0)
  {
    // copy distribution
    _dist = new vector < float >;

    // cerr << "vector<int> input: ";
    for (int i = 0; i < (int) pArray->size(); i++)
    {
      _dist->push_back((float) (*pArray)[i]);
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
  Class:        
  Arguments: 
  
  Description:
  
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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
    cerr <<
      "--Error-- Distribution::TakeSample(): Trying to sample from empty distribution -< return -1."
      << endl;
    return sample;
  }

  // get random number
  float cumsum = 0;
  float value = _rn.GetUniform();
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
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
Distribution::SelectMax()
{
  int m = 0;
  if (_dist == NULL)
  {
    return -1;
  }
  for (int i = 0; i < (int) _dist->size(); i++)
  {
    // cerr << (*_dist)[i] << "|";
    m = ((*_dist)[m] < (*_dist)[i]) ? i : m;
  }
  // cerr << endl;
  return m;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
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

  vector < float >::iterator iter = _dist->begin();
  for (int i = 0; i < length; i++)
    *(iter + start + i) = 0;

  // renormalize distribution
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
bool
Distribution::_Normalize()
{

  // do nothing if distribution is empty
  if (_dist == NULL)
  {
    cerr << "--Error-- Distribution::_Normalize(): _dist is NULL." << endl;
    return false;
  }

  float sum = 0;
  vector < float >::iterator iter = _dist->begin();

  // compute sum over array
  // cerr << "Normalizing: ";
  for (; iter != _dist->end(); iter++)
  {
    sum += *iter;
    // cerr << sum << " ";
  }
  // cerr << endl;

  if (sum == 0)
  {
    cerr << "--Error-- Distribution::_Normalize():  sum is equal to 0." <<
      endl;
    return false;
  }

  // update values
  // cerr << "Normalized values: ";
  for (iter = _dist->begin(); iter != _dist->end(); iter++)
  {
    (*iter) = (*iter) / sum;
    // cerr << (*iter) << " ";
  }
  // cerr << endl;

  return true;
}

#include "MaskVector.h"

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
MaskVector::MaskVector(int l)
{
  // new mask vector
  _mask = new vector < int >(l, 1);
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
MaskVector::~MaskVector()
{
  if (_mask != NULL)
    delete _mask;
  _mask = 0;
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
  MaskVector::ResetMask()
{
  vector < int >::iterator
    iter = _mask->begin();

  while (iter != _mask->end())
  {
    (*iter) = 1;
    iter++;
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
  MaskVector::ClearMask()
{
  vector < int >::iterator
    iter = _mask->begin();

  while (iter != _mask->end())
  {
    (*iter) = 0;
    iter++;
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
  MaskVector::UpdateMask(int start, int w, int value)
{
  int end = start + w;
  if (start < 0)
    start = 0;
  if (start < 0)
    start = 0;
  if (end > (int) _mask->size())
    end = (int) _mask->size();

  for (int i = start; i < end; i++)
    (*_mask)[i] = value;

  //~ for (uint i = 0; i < _mask->size(); i++)
    //~ cerr << (*_mask)[i] << " ";
  //~ cerr << endl;
  
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
  MaskVector::GetValueAt(int i)
{
  if (i >= 0 && i < (int) _mask->size())
  {
    return (*_mask)[i];
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
MaskVector::SetValueAt(int i, int value)
{
  if (i >= 0 && i < (int) _mask->size() && (value == 1 || value == 0))
  {
    (*_mask)[i] = value;
  }
  else
  {
    cerr <<
      "--Warning-- MaskVector::SetValueAt() wrong arguments -> nothing changed."
      << endl;
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
MaskVector::ReverseMask()
{
  int l = (int) _mask->size();
  // create new mask
  MaskVector *pRevMask;
  pRevMask = new MaskVector(l);

  for (int i = 0; i < l; i++)
  {
    pRevMask->SetValueAt(l - i - 1, (*_mask)[i]);
  }

  return pRevMask;
}

#include "InstanceMap.h"
#include "inclusive.h"


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
InstanceMap::InstanceMap()
{
  // create empty list
  _pInstanceList = new list < Instance * >;
  _numberOfInstances = 0;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
InstanceMap::~InstanceMap()
{

  // delete all instances in the list
  list < Instance * >::iterator lIter = _pInstanceList->begin();
  while (lIter != _pInstanceList->end())
  {
    delete(*lIter);
    lIter++;
  }

  // delete the list itself
  delete _pInstanceList;
  _pInstanceList = NULL;
  _numberOfInstances = 0;
}



/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
InstanceMap *
InstanceMap::AddInstance(Instance * pInst)
{

  // add instance at end of list
  _pInstanceList->push_back(pInst);

  // augment number of instances
  _numberOfInstances++;
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
  InstanceMap::ClearMap()
{

  // delete all instances in the list
  list < Instance * >::iterator lIter = _pInstanceList->begin();
  while (lIter != _pInstanceList->end())
  {

    // cerr << "InstanceMap::ClearMap() clear content of iterator" << endl;
    if ((*lIter != NULL))
      delete(*lIter);
    (*lIter) = NULL;
    lIter++;
  }

  // cerr << "InstanceMap::ClearMap() size of list = " << _pInstanceList->size() << endl;

  /* 
   * lIter = _pInstanceList->begin();
   * while ( lIter != _pInstanceList->end() )
   * {          
   * cerr << "InstanceMap::ClearMap() delete Instance from list" << endl; 
   * cerr << "InstanceMap::ClearMap() size of list before erase = " << _pInstanceList->size() << endl;
   * _pInstanceList->erase(lIter);
   * cerr << "InstanceMap::ClearMap() size of list after erase = " << _pInstanceList->size() << endl;
   * lIter++;
   * }  
   */

  // cerr << "InstanceMap::ClearMap() delete current map." << endl;
  delete
    _pInstanceList;

  // cerr << "InstanceMap::ClearMap() create new map." << endl;
  _pInstanceList = new list < Instance * >;
  _numberOfInstances = 0;
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
InstanceMap *
InstanceMap::ReplaceInstance(Instance * oldInst, Instance * newInst)
{

  // locate old instance
  list < Instance * >::iterator lIter = _pInstanceList->begin();
  while (lIter != _pInstanceList->end() && (*lIter) != oldInst)
  {
    lIter++;
  }
  if (lIter != _pInstanceList->end())
  {
    // delete old instance
    delete(*lIter);

    // set new instance
    (*lIter) = newInst;
  }
  else
  {
    // at the end of the list, so add element
    _pInstanceList->push_back(newInst);
    _numberOfInstances++;
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
InstanceMap *
InstanceMap::RemoveInstanceFromSequence(SequenceObject * pSeq)
{
  list < Instance * >::iterator iter = _pInstanceList->begin();
  while (iter != _pInstanceList->end())
  {
    if ((*iter)->ParentSequence() == pSeq)
    {
      delete(*iter);
      _pInstanceList->erase(iter);
    }
    else
    {
      iter++;
    }
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
InstanceMap *
InstanceMap::RemoveInstance(Instance * pInst)
{
  list < Instance * >::iterator iter = _pInstanceList->begin();
  while (iter != _pInstanceList->end())
  {
    if ((*iter) == pInst)
    {
      delete(*iter);
      _pInstanceList->erase(iter);
      break;                    // element is erased leave the while loop
    }
    else
    {
      // next element
      iter++;
    }
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
InstanceMap *
InstanceMap::erase(list < Instance * >::iterator iter)
{
  if ((*iter) != NULL)
    delete(*iter);
  _pInstanceList->erase(iter);
  return this;
}

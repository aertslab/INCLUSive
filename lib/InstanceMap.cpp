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

  if ( _pInstanceList != NULL )
  {    
    // delete all instances in the list
    list < Instance * >::iterator lIter = _pInstanceList->begin();
    while (lIter != _pInstanceList->end())
    {
      if ( *lIter != NULL )
        delete(*lIter);
      lIter++;
    }

    // cerr << "InstanceMap::ClearMap() delete current map." << endl;
    delete
      _pInstanceList;
  }
  
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


InstanceMap* 
InstanceMap::ExtendLeft(int pos){
  list<Instance *>::iterator lIter = _pInstanceList->begin();
  while (lIter != _pInstanceList->end())
  {
    if ( !(*lIter)->CheckLeftExtension(pos) )
      return this; // return without adjusting anything
    lIter++;
  }
  
  // now the instances can be extended
  lIter = _pInstanceList->begin();
  while (lIter != _pInstanceList->end())
  {
    (*lIter)->ExtendLeft(pos);
    lIter++;
  }
  
  return this;
} 


InstanceMap* 
InstanceMap::ExtendRight(int pos){
  list < Instance * >::iterator lIter = _pInstanceList->begin();
  while (lIter != _pInstanceList->end())
  {
    if ( !(*lIter)->CheckRightExtension(pos) )
      return this; // return without adjusting anything
    lIter++;
  }
  
  // now the instances can be extended
  lIter = _pInstanceList->begin();
    while (lIter != _pInstanceList->end())
  {
    (*lIter)->ExtendRight(pos);
    lIter++;
  }
  
  return this;
}

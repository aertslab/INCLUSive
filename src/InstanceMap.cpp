#include "InstanceMap.h"
#include "inclusive.h"


/******************************************************************************
  Method:       new
  Class:        InstanceMap 
  Arguments:    none
  
  Description:  constructor
	              creates an empty list of Instance objects
  
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
  Method:       delete
  Class:        InstanceMap
  Arguments:    none
  
  Description:  destructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
InstanceMap::~InstanceMap()
{

  // first delete all individual instances in the list
  list < Instance * >::iterator lIter = _pInstanceList->begin();
  while (lIter != _pInstanceList->end())
  {
	if ( *lIter != NULL )
      delete(*lIter);
    lIter++;
  }

  // delete the list itself
  delete _pInstanceList;
  _pInstanceList = NULL;
  _numberOfInstances = 0;
}



/******************************************************************************
  Method:       AddInstance
  Class:        InstanceMap
  Arguments:    Instance * pInst
  
  Description:  add a new instance to the list
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
InstanceMap *
InstanceMap::AddInstance(Instance * pInst)
{
  // add instance at end of list
	if ( pInst != NULL )
	{
		_pInstanceList->push_back(pInst);

		// augment number of instances
		_numberOfInstances++;
	}
  return this;
}



/******************************************************************************
  Method:       ClearMap
  Class:        InstanceMap
  Arguments:    none
  
  Description:  clear all instances in the list and 
                start with a new empty list
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
  InstanceMap::ClearMap()
{

	// first delete the current list
  if ( _pInstanceList != NULL )
  {    
    // delete all individual instances in the list
    list < Instance * >::iterator lIter = _pInstanceList->begin();
    while (lIter != _pInstanceList->end())
    {
      if ( *lIter != NULL )
        delete (*lIter);
      lIter++;
    }
    delete _pInstanceList;
  }
  
	// create a new empty list
  _pInstanceList = new list < Instance * >;
  _numberOfInstances = 0;
  return;
}



/******************************************************************************
  Method:       ReplaceInstance
  Class:        InstanceMap
  Arguments:    Instance * oldInst, Instance * newInst
  
  Description:  replace one instance in the list with another instance
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
InstanceMap *
InstanceMap::ReplaceInstance(Instance * oldInst, Instance * newInst)
{

  // locate old instance
  //~ list < Instance * >::iterator lIter = _pInstanceList->begin();
  //~ while (lIter != _pInstanceList->end() && (*lIter) != oldInst)
  //~ {
    //~ lIter++;
  //~ }
	list < Instance * >::iterator lIter = find(_pInstanceList->begin(),_pInstanceList->end(),oldInst);

  if (lIter != _pInstanceList->end())
  {
    // delete old instance
		if ( (*lIter) != NULL )
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
  Method:       RemoveInstanceFromSequence
  Class:        InstanceMap
  Arguments:    SequenceObject * pSeq
  
  Description:  remove instances from sequence pSeq from the list
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
InstanceMap *
InstanceMap::RemoveInstanceFromSequence(SequenceObject * pSeq)
{
	// iterate over all instances
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
  Method:       RemoveInstance
  Class:        InstanceMap
  Arguments:    Instance * pInst
  
  Description:  remove a specific instance from the list
  
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
  Method:       Erase
  Class:        InstanceMap
  Arguments:    list<Instance *>::iterator iter
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
InstanceMap *
InstanceMap::Erase(list<Instance *>::iterator iter)
{
  if ((*iter) != NULL)
    delete(*iter);
  _pInstanceList->erase(iter);
  return this;
}


/******************************************************************************
  Method:       ExtendLeft
  Class:        InstanceMap
  Arguments:    int p
  
  Description:  try to extend all instances in the list p positions
                to the left
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
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


/******************************************************************************
  Method:       ExtendRight
  Class:        InstanceMap
  Arguments:    int pos
  
  Description:  try to extend all instances in the list p positions
                to the right
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
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

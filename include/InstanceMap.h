#ifndef instancemap_include_declared
#define instancemap_include_declared

#include "Instance.h"
#include <list>

class InstanceMap {
 private:
  int _numberOfInstances;
  int _iterCounter;
  list<Instance *> *_pInstanceList;
  
 public:
  // constructor
  InstanceMap();
  
  // destructor
  ~InstanceMap();
  
  // inspectors
  int NumberOfInstances(){return _numberOfInstances;};
  list<Instance *>::iterator begin(){return _pInstanceList->begin();};
  list<Instance *>::iterator end(){return _pInstanceList->end();};

  // adaptors
  InstanceMap * AddInstance(Instance *pInst);
  InstanceMap * RemoveInstance(Instance *pInst);
  InstanceMap * ReplaceInstance(Instance *oldInst, Instance *newInst);
  InstanceMap * RemoveInstanceFromSequence(SequenceObject *pSeq);
  InstanceMap * erase(list<Instance *>::iterator iter);
  void ClearMap();
  
};

#endif

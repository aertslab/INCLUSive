#include "Instance.h"

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Instance::Instance(SequenceObject * pSeqObj, strand_modes s, int start,
                   int length)
{
  // set parameters
  if (start < 0 || length < 0)
  {
    _start = 0;
    _length = 0;
    _score = 0;
    _pSite = NULL;
    _pSiteString = NULL;
    _pSeqName = NULL;
    return;
  }

  // set internal parameters
  _strand = s;
  _start = start;
  _length = length;
  _pParentSequence = pSeqObj;
  _pSeqName = new string(*(pSeqObj->GetID()));

  // initial score is 0
  _score = 0;

  // get site
  _pSite = _pParentSequence->GetSubSequenceVector(s, start, length);

  // translate sequence vector to 
  if (_pSite == NULL)
  {
    _pSiteString = NULL;
    cerr <<
      "-Error-- Instance(): Unable to get site from sequence => No instance created "
      << endl;
  }
  else
  {

    _pSiteString = new string;

    for (uint i = 0; i < _pSite->size(); i++)
    {
      switch ((*_pSite)[i])
      {
      case 0:
        _pSiteString->append("A");
        break;
      case 1:
        _pSiteString->append("C");
        break;
      case 2:
        _pSiteString->append("G");
        break;
      case 3:
        _pSiteString->append("T");
        break;
      default:
        _pSiteString->append("N");
      }
    }
    // cerr << "Instance " << index << " => created string: " << *_pSiteString << endl;
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
Instance::Instance(Instance& copyInstance)
{
    // set internal parameters
  _strand = copyInstance.Strand();
  _start = copyInstance.Start();
  _length = copyInstance.Length();
  _pParentSequence = copyInstance.ParentSequence();
  _pSeqName = new string(*(copyInstance.ParentSequence()->GetID()));
  _score = copyInstance.Score();
  _pSite = _pParentSequence->GetSubSequenceVector(_strand, _start, _length);
  _pSiteString = new string(*(copyInstance.PrintSite()));
  
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Instance::~Instance()
{
  // reset pointer to string (should be defined outside this scope);
  // _pId = NULL;

  if (_pSite != NULL)
    delete _pSite;
  _pSite = NULL;

  if (_pSiteString != NULL)
    delete _pSiteString;
  _pSiteString = NULL;

  if ( _pSeqName != NULL )
    delete _pSeqName;
  _pSeqName = NULL;
  
}

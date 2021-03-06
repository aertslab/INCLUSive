#include "Instance.h"

/******************************************************************************
  Method:       new 
  Class:        Instance
  Arguments:    SequenceObject * pSeqObj, strand_modes s, int start, int length
  
  Description:  constructor
							  create new instance from the given sequence
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
  Method:       new
  Class:        Instance
  Arguments:    Instance& copyInstance
  
  Description:  copy constructor
  
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
  Class:        Instance
  Arguments:    none
  
  Description:  destructor
  
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



/******************************************************************************
  Method:       CheckLeftExtension
  Class:        Instance
  Arguments:    int p
  
  Description:  Check if the instance is extensible p positions to the left
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
bool Instance::CheckLeftExtension(int pos){
  if ( _start - pos < 0 ){
    cerr << "--Warning-- Instance::ExtendLeft(): Unable to extend site to left." << endl;
    return false;
  }
  
  return true;
}


/******************************************************************************
  Method:       ExtendLeft
  Class:        Instance
  Arguments:    int p
  
  Description:  Extend the instance p positions to the left
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Instance* Instance::ExtendLeft(int pos){
  // change start position ans site length
  _start = _start - pos;
  _length += pos;
  
  // update site
  if (_pSite != NULL)
    delete _pSite;
  // cerr << "DEBUG: Instance::ExtendLeft(int pos) " << pos << " -> " << _start << " " << _length << endl;
  _pSite = _pParentSequence->GetSubSequenceVector(_strand, _start, _length);
  
  // update site string
  _pSiteString->resize(_length);
  for (int i = 0; i < _length; i++)
    {
      switch ((*_pSite)[i])
      {
      case 0:
        _pSiteString->replace(i,1,"A");
        break;
      case 1:
        _pSiteString->replace(i,1,"C");
        break;
      case 2:
        _pSiteString->replace(i,1,"G");
        break;
      case 3:
        _pSiteString->replace(i,1,"T");
        break;
      default:
        _pSiteString->replace(i,1,"N");
      }
    }
  
  return this;
}


/******************************************************************************
  Method:       CheckRightExtension
  Class:        Instance
  Arguments:    int p
  
  Description:  Check if the instance is extensible p positions to the right
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
bool Instance::CheckRightExtension(int pos){
  if ( _start + _length + pos > _pParentSequence->Length() ){
    cerr << "--Warning-- Instance::ExtendRight(): Unable to extend site to right." << endl;
    return false;
  }
  
  return true;
}

  
/******************************************************************************
  Method:       ExtendRight
  Class:        Instance
  Arguments:    int p
  
  Description:  Extend the instance p positions to the right
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
Instance*  Instance::ExtendRight(int pos){
  // change length of the site
  _length += pos;
  
  // update site
  if (_pSite != NULL)
    delete _pSite;

  // cerr << "DEBUG: Instance::ExtendRight(int pos) " << pos << " -> " << _start << " " << _length << endl;
  _pSite = _pParentSequence->GetSubSequenceVector(_strand, _start, _length);

  for (int i = _length-pos; i < _length; i++)
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
  return this;
}

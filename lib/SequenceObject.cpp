#include "SequenceObject.h"
#include "utilities.h"

// constructor
SequenceObject::SequenceObject(const string * pSeqStr)
{

  _length = pSeqStr->size();

  // set sequence vector
  // define new sequence vector
  _sequenceVector = new vector < int >;
  _CreateSequenceVector(pSeqStr);
  _revSequenceVector = new vector < int >(_length);
  _SetReverseComplement();

  // set new identifiers
  _id = "unknown";
  _description = "unknown";

}


// constructor
SequenceObject::SequenceObject(Sequence * pSeqVector)
{

  _length = pSeqVector->size();

  // set sequence vector
  // define new sequence vector
  _sequenceVector = new vector < int >;
  vector < int >::iterator sIter = pSeqVector->begin();
  for (; sIter != pSeqVector->end(); sIter++)
  {
    _sequenceVector->push_back(*(sIter));
  }

  // set new identifiers
  _id = "unknown";
  _description = "unknown";

  // create reverse complement
  _revSequenceVector = new vector < int >(_length);
  _SetReverseComplement();

}


// destructor
SequenceObject::~SequenceObject()
{
  // free sequence vectors
  if (_sequenceVector != NULL)
    delete _sequenceVector;
  _sequenceVector = NULL;

  if (_revSequenceVector != NULL)
    delete _revSequenceVector;
  _revSequenceVector = NULL;

}


// inspectors
int
SequenceObject::Length()
{
  return _length;
}

string *
SequenceObject::GetID()
{
  return &_id;
}

string *
SequenceObject::GetDescription()
{
  return &_description;
}

string *
SequenceObject::GetSequenceString(strand_modes s)
{
  return _TranslateSequenceVector(s, 0, _length);
}


int
SequenceObject::GetNucleotideAt(strand_modes s, int pos)
{
  if (pos < 0 || pos >= _length)
  {
    return -1;
  }
  else if (s == plus_strand)
  {
    return (*_sequenceVector)[pos];
  }
  else if (s == minus_strand)
  {
    return (*_revSequenceVector)[pos];
  }
  else
  {
    return -1;
  }
}

Sequence *
SequenceObject::GetSubSequenceVector(strand_modes s, int start, int length)
{
  vector < int >::iterator iterSeq;
  int ll = _sequenceVector->size();

  if (start < 0 || ll < start + length)
  {
    cerr << "Warning: start not within range of the sequence." << endl;
    return NULL;
  }

  // define a vector to store the results
  vector < int >*pRes = new vector < int >(length);

  if (s == 1)
  {
    iterSeq = _sequenceVector->begin();
  }
  else
  {
    iterSeq = _revSequenceVector->begin();
  }

  for (int i = 0; i < length; i++)
    (*pRes)[i] = *(iterSeq + start + i);

  return pRes;

}


string *
SequenceObject::GetSubSequenceString(strand_modes s, int start, int length)
{
  return _TranslateSequenceVector(s, start, length);
}


bool
SequenceObject::CheckSequence()
{
  bool bCheck = true;

  vector < int >::iterator iterSeq = _sequenceVector->begin();

  for (; iterSeq != _sequenceVector->end(); iterSeq++)
  {
    if (*iterSeq == -1)
    {
      bCheck = false;
      break;
    }
  }
  return bCheck;

}


// adaptors
SequenceObject *
SequenceObject::SetID(string * idname)
{
  _id = *idname;
  return this;
}


SequenceObject *
SequenceObject::SetDescription(string * desc)
{
  _ParseDescriptionHeader(desc);
  return this;
}


string *
SequenceObject::_TranslateSequenceVector(strand_modes s, int start,
                                         int length)
{
  vector < int >::iterator iterB;
  // check parameters

  // create new string
  string *pStr = new string("");

  if (s == 1)
  {
    iterB = _sequenceVector->begin();
  }
  else
  {
    iterB = _revSequenceVector->begin();
  }

  for (int i = start; i < start + length; i++)
  {
    if (i >= length)
      break;

    switch (*(iterB + i))
    {
    case 0:
      pStr->append("A");
      break;
    case 1:
      pStr->append("C");
      break;
    case 2:
      pStr->append("G");
      break;
    case 3:
      pStr->append("T");
      break;
    default:
      pStr->append("N");
    }
  }
  return pStr;
}


void
SequenceObject::_CreateSequenceVector(const string * pSeqStr)
{
  int i = 0;

  // transform the sequence (A=0,C=2,G=3,T=4)
  while ( i < (int)pSeqStr->length())
  {
    switch (pSeqStr->at(i))
    {
    case 'a':
    case 'A':
      _sequenceVector->push_back(0);
      break;
    case 'c':
    case 'C':
      _sequenceVector->push_back(1);
      break;
    case 'g':
    case 'G':
      _sequenceVector->push_back(2);
      break;
    case 't':
    case 'T':
      _sequenceVector->push_back(3);
      break;
    default:
      _sequenceVector->push_back(-1);
    }
    //cerr << pSeqStr->at(i) << (* _sequenceVector)[i];
    i++;
  }

  //cerr << endl;
  return;
}

void
SequenceObject::_SetReverseComplement()
{
  int pos = 0;
  for (int i = 0; i < _length; i++)
  {
    pos = _length - 1 - i;
    switch ((*_sequenceVector)[i])
    {
    case 0:
      (*_revSequenceVector)[pos] = 3;
      break;
    case 1:
      (*_revSequenceVector)[pos] = 2;
      break;
    case 2:
      (*_revSequenceVector)[pos] = 1;
      break;
    case 3:
      (*_revSequenceVector)[pos] = 0;
      break;
    default:
      (*_revSequenceVector)[pos] = -1;
    }
  }

  return;
}


vector < int >*
_TranslateSequenceString(const string * seqStr)
{

  vector < int >*pRes;
  pRes = new vector < int >;

  int i = 0;

  // transform the sequence (A=0,C=2,G=3,T=4)
  while (i < (int)seqStr->length())
  {
    switch (seqStr->at(i))
    {
    case 'a':
    case 'A':
      pRes->push_back(0);
      break;
    case 'c':
    case 'C':
      pRes->push_back(1);
      break;
    case 'g':
    case 'G':
      pRes->push_back(2);
      break;
    case 't':
    case 'T':
      pRes->push_back(3);
      break;
    default:
      pRes->push_back(-1);
    }
    i++;
  }

  return pRes;
}


void
SequenceObject::_ParseDescriptionHeader(string * descStr)
{
  string tokens(" \t");
  string::size_type pos = 0;

  if ((pos = descStr->find_first_of(tokens)) == string::npos)
  {
    _id = *descStr;
    _description = *descStr;
  }
  else
  {
    // int l = descStr.size();
    _id = descStr->substr(1, pos-1);
    _description = *descStr;
  }

  // remove trailing white space from description
  tokens = "> \t";
  while ((pos = _id.find_first_of(tokens)) == 0)
  {
    _id.erase(pos, 1);
  }


  tokens = "\r\f\n";
  while ((pos = _id.find_first_of(tokens)) != string::npos)
  {
    _id.erase(pos, 1);
  }

  while ((pos = _description.find_first_of(tokens)) != string::npos)
  {
    _description.erase(pos, 1);
  }

  return;
}

Sequence::iterator 
SequenceObject::begin(strand_modes s)
{
  if (s == plus_strand)
  {
    return _sequenceVector->begin();
  }
  else
  {
    return _revSequenceVector->begin();
  }
}


Sequence::iterator SequenceObject::end(strand_modes s)
{
  if (s == plus_strand)
  {
    return _sequenceVector->end();
  }
  else
  {
    return _revSequenceVector->end();
  }
}

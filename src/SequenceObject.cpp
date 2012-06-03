#include "SequenceObject.h"
#include "utilities.h"

/****************************************************************************
  Method:       new
  Class:        SequenceObject
  Arguments:    const string *pSeqStr
  
  Description:  constructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       new
  Class:        SequenceObject
  Arguments:    Sequence * pSeqVector
  
  Description:  constructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       delete
  Class:        SequenceObject
  Arguments:    none
  
  Description:  destructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       GetSequenceString
  Class:        SequenceObject
  Arguments:    strand_modes s
  
  Description:  returns a string that contains the sequence
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
string *
SequenceObject::GetSequenceString(strand_modes s)
{
  return _TranslateSequenceVector(s, 0, _length);
}


/****************************************************************************
  Method:       GetNucleotideAt
  Class:        SequenceObject
  Arguments:    strand_modes s, int pos
  
  Description:  gives the nucleotide at position pos of strand s
                values are (A=0,C=1,G=2,T=3, or -1)
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       GetSubSequenceVector
  Class:        SequenceObject
  Arguments:    strand_modes s, int start, int length
  
  Description:  returns the pointer to a vector that contains
                the subsequence from start to start+length-1 
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
Sequence *
SequenceObject::GetSubSequenceVector(strand_modes s, int start, int length)
{
  vector < int >::iterator iterSeq;
  int ll = _sequenceVector->size();

  if (start < 0 || ll < start + length)
  {
    cerr << "Warning: start not within range of the sequence " << _id << ": " << ll << " = " << start << "<->" <<start + length << endl;
    return NULL;
  }

  // define a vector to store the results
  vector<int>* pRes = new vector<int>(length);

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


/****************************************************************************
  Method:       GetSubSequenceString
  Class:        SequenceObject
  Arguments:    strand_modes s, int start, int length
  
  Description:  returns the pointer to a string that contains
                the subsequence from start to start+length-1 
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
string *
SequenceObject::GetSubSequenceString(strand_modes s, int start, int length)
{
  return _TranslateSequenceVector(s, start, length);
}


/****************************************************************************
  Method:       CheckSequence
  Class:        SequenceObject
  Arguments:    none

  Description:  check if the sequence has the correct symbols
	              returns -1 if one none ACGT symbol is found

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       _TranslateSequenceVector
  Class:        SequenceObject
  Arguments:    strand_modes s, int start, int length
  
  Description:  internal method
                returns the pointer to a string that contains
                the translated subsequence from start to start+length-1 
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       _CreateSequenceVector
  Class:        SequenceObject
  Arguments:    const string * pSeqStr
  
  Description:  internal method to translate a sequence vector in 
                numerical format

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       _SetReverseComplement
  Class:        SequenceObject
  Arguments:    none
  
  Description:  internal method
								to create the reverse complement of the sequence vector

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       _TranslateSequenceString
  Class:        SequenceObject
  Arguments:    strand_modes s, int start, int length
  
  Description:  internal method
                returns the pointer to a vector that contains
                the translated subsequence from start to start+length-1 
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
vector<int>*
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


/****************************************************************************
  Method:       _ParseDescriptionHeader
  Class:        SequenceObject
  Arguments:    strand_modes s, int start, int length
  
  Description:  internal method to parse header in fasta file
                format is '>id description'
                everything between > and first white space is id
                rest is descriptor

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       begin
  Class:        SequenceObject
  Arguments:    strand_modes s
  
  Description:  begin iterator 

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
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


/****************************************************************************
  Method:       end
  Class:        SequenceObject
  Arguments:    strand_modes s
  
  Description:  end iterator 

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
****************************************************************************/
Sequence::iterator 
SequenceObject::end(strand_modes s)
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

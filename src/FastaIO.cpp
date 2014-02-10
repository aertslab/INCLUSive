// 21 july 2009 - 3.1.5

#include "FastaIO.h"

#include <iostream>


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26 
                rev 2011/07/03 Marleen Claeys : 
                    allow phylofasta file as input in MS. 
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
FastaIO::FastaIO(string * fileName)
{
  // open file for reading
  _ifs.open(fileName->c_str(), ios::in);

  if (!_ifs)
  {
    // file is not open for reading
    cerr << "--Error-- FastaIO::FastaIO() unable to open file: "
      << *fileName << " for reading!" << endl;
    _isOpen = false;
    _hasNext = false;
    /* 
     * the program will not stop here
     * the user should check with IsOpen() in main program
     */
  }
  else
  {
    _isOpen = true;

    // read lines untill first line starting with a '>' is found 
    do
    {
      getline(_ifs, _pLine, '\n');
      // cerr << "DEBUG: " << _pLine << endl;
    }
    //while (_pLine[0] != '>' && !_ifs.eof()); // but skip lines with double >> (refers to gene-group)
    while (!(_pLine[0] == '>' && _pLine[1] != '>') && !_ifs.eof());

    if (!_ifs.eof())
    {
      _hasNext = true;
      // cerr << "DEBUG: _hasNext = true " << endl;
    }
    else
    {
      _hasNext = false;
      // cerr << "DEBUG: _hasNext = false " << endl;
    }
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
FastaIO::~FastaIO()
{
  // close file
  if (_ifs.good())
    _ifs.close();

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
FastaIO::Close()
{
  _ifs.close();
  _isOpen = false;
  return;
}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  Date:         MC_2011/07/03 allow phylofasta file as input in MS.
  
******************************************************************************/
SequenceObject *
FastaIO::NextSequence()
{

  if (_isOpen && _hasNext)
  {

    // create empty string to store sequence
    string *pSeqStr = new string("");

    string *pDescStr = new string(_pLine);
    
	  //cerr << "DEBUG:FastaIO::NextSequence:: " << *pDescStr;

    while (!_ifs.eof())
    {
      // read next line
      getline(_ifs, _pLine, '\n');

      // move to next line if line is empty and skip >>gene-group lines 
      if (_pLine.empty() || (_pLine[0] == '>' && _pLine[1] == '>'))
        continue;

      if (_pLine[0] == '>')
      {
        break;
      }
      else
      {
        string::size_type pos = 0;
        string tokens = " \t\n\r\f";
        while ((pos = _pLine.find_first_of(tokens)) != string::npos)
        {
          _pLine.erase(pos, 1);
        }
        // append line to seqStr;
        pSeqStr->append(_pLine);

      }
    }

    if (!_ifs.eof())
    {
      _hasNext = true;
    }
    else
    {
      _hasNext = false;
    }

    // create new sequence object
    // cerr << *seqStr << endl;

    SequenceObject *pSeqObj = new SequenceObject(pSeqStr);
    pSeqObj->SetDescription(pDescStr);

    delete pSeqStr;
    delete pDescStr;

    return pSeqObj;

  }
  else
  {
    // sequence was not open anymore 
    return NULL;
  }
}


/******************************************************************************
  Description:  read line of format >value_value_value_value
  Author_Date:  MC_2009/07/21
******************************************************************************/
double *
FastaIO::ReadDirichlet()
{
  // only this line is read, cursor is not moved to next line
    
  double * output = new double[4];
  string::size_type pos;
  string * input = new string(_pLine);
  input->append("_");
  double v;
  input->erase(0,1); // erase >
  for(int i = 0; i < 4; i++)
  {
    pos = input->find_first_not_of("0.123456789");
      if ( pos == string::npos )         
    { cerr << "Error--FastaIO:ReadDirichlet(): format not correct.\n";
      delete[] output; output = NULL; break; 
    }
    // convert to double
    istringstream istr(input->substr(0,pos));
    istr >> v; 
    if ( v < 0)         
    { cerr << "ERROR--FastaIO:ReadDirichlet() : values should be > 0.\n"; 
      delete[] output; output = NULL; break; 
    }
    output[i] = v;
    input->erase(0,pos+1); 
  }
  delete input; input = NULL;

  return output;
}

/******************************************************************************
  Description:  return id in line of format >id
                and put pointer on the next line
  Author_Date:  MC_2012/09/19
******************************************************************************/
string *
FastaIO::ReadSeqID()
{
  // this line is read, cursor is moved to next line
  string * readid = NULL;

  if (_isOpen && _hasNext)
  {
    // store _pLine into string
    readid = new string(_pLine);
    readid->erase(0,1); // erase >
    if ( readid->length() < 1) 
    { delete readid; return NULL;} // no id was provided

    // read next line
    _pLine = ""; // reset
    if (!_ifs.eof())
    {
      getline(_ifs, _pLine, '\n'); // this line should contain probs-data
    }
    _hasNext = false; // _pLine does not start with > now
  }
  return readid;
}
/******************************************************************************
  Description:  read all the probs (L, or 2L) available in _pLine
                the values are separated by a space
                ! the last L-w values are (if not yet so) overwritten to zero !
                and put pointer on the next line (>)
				Rescale all entries to value range [0.1-0.9]
  Author_Date:  MC_2014/01/12
******************************************************************************/
ScoreVector *
FastaIO::LoadPspData(int seqL, int wL, bool skip)
{
  ScoreVector * pData = NULL;
if (!skip) // else skip processing this line and put ready on next >line
{
  // _pLine format should be e.g. 0.4 0.6 0.8 ... 0.0 0.0 
  string::size_type pos;
  pos = _pLine.find_first_not_of(" \t0.123456789\n\r\f");
  if ( pos != string::npos )
  { return NULL;}
  // load the probs into Scorevector
  int i = 0; // count how many entries loaded
  double Max = -1; double Min = 100000; // track to enable rescaling [0.1-0.9]
  double prob;
  pData = new ScoreVector;
  // now load data
  while (!_pLine.empty())
  { 
    pos = _pLine.find_first_not_of("0.123456789");
    istringstream istr(_pLine.substr(0,pos)); // convert to double
    istr >> prob;
    _pLine.erase(0,pos); // erase prob
    while ( _pLine.find_first_of("\t ") == 0) _pLine.erase(0,1); // erase space
    if ((i > seqL-wL && i < seqL) || (i > 2*seqL-wL && i < 2*seqL))
    { pData->push_back(0);} // overwrite to zero (masked positions at seq-end)
    else 
    { pData->push_back(prob);
      if (prob > Max){Max = prob;}
      if (prob < Min){Min = prob;}
    }
    i++;
  }
  // check if total pData equals L or 2L
  if ( i != seqL && i != 2*seqL)
  { cerr << "--Error--PSP: nbr of entries = " << i << " while L/2L = " << seqL << "/" << 2*seqL << endl;
    delete pData; pData = NULL; return NULL;}
  // rescale all valid loaded entries between [0.1- 0.9]
  if (Max != Min) // else all entries are the same (uniform psp)
  {
    for (int j = 0; j < i; j++)
    { 
      if (!((j > seqL-wL && j < seqL) || (j > 2*seqL-wL && j < 2*seqL)))
      {(*pData)[j] = 0.1 + 0.8*((*pData)[j]-Min)/(Max-Min);}
    }
  }
  else // set all entries to '1' (uniform psp)
  {
    for (int j = 0; j < i; j++)
    { if (!((j > seqL-wL && j < seqL) || (j > 2*seqL-wL && j < 2*seqL)))
      {(*pData)[j] = 1;}
    }
  }
}


  // set pointer ready on next line
  _pLine = "";
  while (!_ifs.eof())
  {
    // read next line
    getline(_ifs, _pLine, '\n');
    // move to next line if line is empty and skip >>gene-group lines 
    if (_pLine.empty() || (_pLine[0] == '>' && _pLine[1] == '>'))
      continue;
    if (_pLine[0] == '>')
    {
      break;
    }
  }
  if (!_ifs.eof())
  {
    _hasNext = true;
  }
  else
  {
    _hasNext = false;
  }

  return pData;
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
FastaIO::IsOpen()
{
  return _isOpen;
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
FastaIO::HasNext()
{
  return _hasNext;
}

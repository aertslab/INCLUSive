#include "FastaIO.h"

#include <iostream>


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
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
    while (_pLine[0] != '>' && !_ifs.eof());

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
  
******************************************************************************/
SequenceObject *
FastaIO::NextSequence()
{

  if (_isOpen && _hasNext)
  {

    // create empty string to store sequence
    string *pSeqStr = new string("");

    string *pDescStr = new string(_pLine);
    
    // cerr << "DEBUG: " << *pDescStr;

    while (!_ifs.eof())
    {
      // read next line
      getline(_ifs, _pLine, '\n');

      // move to next line if line is empty
      if (_pLine.empty())
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

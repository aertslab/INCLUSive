#include "PWMIO.h"

#include <cstring>
#include <stdio.h>

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWMIO::PWMIO(string * fileName, int type)
{
  // define string to hold data
  _pLine = new string;

  if (type == READ)
  {

    // set writing parameters off
    _isOpenWriting = false;

    // open file for reading 
    _ifs.open(fileName->c_str(), ios::in);
    if (_ifs)
    {
      // read first line
      getline(_ifs, *_pLine, '\n');
      cerr << "DEBUGLOG: " << *_pLine << endl;
      string is = "#INCLUSive";

      if (_pLine->compare(0, 10, is) == 0)
      {
        // this is a correct matrix file to start with
        _isOpenReading = true;
      }
      else
      {
        // something is wrong with this matrix file
        cerr << "--Error-- PMWIO(): Uncorrect format of matrix file: " <<
          fileName << endl;
        _isOpenReading = false;
        _ifs.close();
      }
    }
    else
    {
      cerr << "--Error-- PMWIO(): Unable to open file: " << fileName <<
        " for reading." << endl;
      _isOpenReading = false;
    }
  }
  else if (type == WRITE)
  {

    // set reading parameters off
    _isOpenReading = false;

    // open file for writing
    _ofs.open(fileName->c_str(), ios::out);
    if (_ofs)
    {
      _isOpenWriting = true;

      // write header line
      _ofs << "#INCLUSive Motif Model v1.0" << endl << "#" << endl;

    }
    else
    {
      cerr << "--Error-- PMWIO(): Unable to open file " << fileName <<
        " for writing." << endl;
      _isOpenWriting = false;
    }
  }
  else
  {
    cerr << "--Error-- PMWIO(): Uncorrect type to open file " << fileName <<
      " Use READ=0 || WRITE=1" << endl;
    _isOpenWriting = false;
    _isOpenReading = false;
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
PWMIO::~PWMIO()
{
  if (_isOpenReading && _ifs)
    _ifs.close();

  if (_isOpenWriting && _ofs)
    _ofs.close();

  if (_pLine != NULL)
    delete _pLine;
  _pLine = NULL;

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
  PWMIO::Close()
{
  if (_isOpenReading && _ifs)
    _ifs.close();

  if (_isOpenWriting && _ofs)
    _ofs.close();

  _isOpenWriting = false;
  _isOpenReading = false;

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
PWM *
PWMIO::ReadMatrix()
{
  int W = 0;
  float fA = 0.25, fC = 0.25, fG = 0.25, fT = 0.25;
  char *pcTag, *pcValue;
  string *sId, *sCons;
  double score = 0;

  // define strings
  pcTag = new char[64];
  pcValue = new char[129];

  sId = new string("");
  sCons = new string("");

  if (!_isOpenReading || !_ifs)
  {
    cerr <<
      "--Error-) PWMIO::ReadMatrix(): Trying to read from file which is not open for reading\n";
    return NULL;
  }

  if (_ifs.eof())
  {
    _isOpenReading = 0;
    delete[]pcTag;
    delete[]pcValue;
    delete sCons;
    delete sId;
    return NULL;
  }

  while (_pLine->empty() || _pLine->length() <= 1 || (*_pLine)[0] == '#')
  {
    int nc = sscanf(_pLine->c_str(), "%s = %s", pcTag, pcValue);

    if (nc == 2)
    {
      if (strcmp(pcTag, "#ID") == 0)
      {
        sId->append(pcValue);
      }
      else if (strcmp(pcTag, "#Consensus") == 0)
      {
        sCons->append(pcValue);
      }
      else if (strcmp(pcTag, "#W") == 0)
      {
        W = atoi(pcValue);
      }
      else if (strcmp(pcTag, "#Score") == 0)
      {
        score = (double) atof(pcValue);
      }
    }
    if (!_ifs.eof())
    {
      getline(_ifs, *_pLine, '\n');
    }
    else
    {
      _isOpenReading = 0;
      delete sCons;
      delete sId;
      delete[]pcTag;
      delete[]pcValue;
      return NULL;
    }
  }

  // clean up temporary variables
  delete[]pcTag;
  delete[]pcValue;


  if (W <= 0)
  {
    cerr <<
      "--Error-- PWMIO::ReadMatrix(): Found non-positive motif length => Returning NULL."
      << endl;
    //clear up some 
    delete sCons;
    delete sId;
    return NULL;
  }

  Matrix pM = new double[W][4];

  for (int i = 0; i < W; i++)
  {
    // check the format of the input line
    int nc = sscanf(_pLine->c_str(), "%f%f%f%f", &fA, &fC, &fG, &fT);
    if (nc == 4)
    {
      pM[i][0] = (double) fA;
      pM[i][1] = (double) fC;
      pM[i][2] = (double) fG;
      pM[i][3] = (double) fT;
    }
    else
    {
      cerr << "--Warning-- PWMIO::ReadMatrix(): line is not as expected: " <<
        *_pLine << endl;
    }

    if (!_ifs.eof())
    {
      getline(_ifs, *_pLine, '\n');
    }
    else
    {
      _isOpenReading = 0;
    }
  }

  // define default snf
  double *snf = new double[4];
  snf[0] = 0.0001;
  snf[1] = 0.0001;
  snf[2] = 0.0001;
  snf[3] = 0.0001;

  // create PWM object
  PWM *pModel = new PWM(W, pM, snf);
  pModel->SetID(sId);
  pModel->SetConsensus(sCons);
  pModel->SetScore(score);

  delete sId;
  delete sCons;
  delete[]pM;
  delete[]snf;
  
  return pModel;

}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
  PWMIO::WriteMatrix(PWM * pMatrix)
{
  int
    W,
    i,
    j;
  string *
  pId, *
    pConsensus;

  if (_ofs && _isOpenWriting)
  {
    // write motif information
    W = pMatrix->Length();
    pId = pMatrix->GetID();
    pConsensus = pMatrix->GetConsensus();

    if (pId == NULL)
    {
      _ofs << "#ID = box" << endl;
    }
    else
    {
      _ofs << "#ID = " << *pId << endl;
    }
    _ofs << "#Score = " << pMatrix->Score() << endl;
    _ofs << "#W = " << W << endl;
    if (pConsensus == NULL)
    {
      _ofs << "#Consensus = not_given" << endl;
    }
    else
    {
      _ofs << "#Consensus = " << *pConsensus << endl;
    }

    // write matrix
    for (i = 0; i < W; i++)
    {
      for (j = 0; j < 4; j++)
      {
        _ofs << pMatrix->GetValueAt(i, j) << "\t";
      }
      _ofs << endl;
    }
    _ofs << endl;

    return 1;

  }
  else
  {
    cerr << "--Error-- PWMIO::WriteMatrix(): File not open for writing" <<
      endl;
    return -1;
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
PWMIO::IsOpen()
{
  if (_isOpenReading || _isOpenWriting)
  {
    return true;
  }
  else
  {
    return false;
  }
}

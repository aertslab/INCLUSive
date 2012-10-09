// 21 july 2009 - 3.1.5

#include "PWMIO.h"

#include <cstring>
#include <stdio.h>

/******************************************************************************
  Method:       new
  Class:        PWMIO
  Arguments:    string * fileName, int type
  
  Description:  opens a file for reading or writing matrices 
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWMIO::PWMIO(string * fileName, int type)
{
  _pError = NULL;
  ostringstream cerrstr;

  string *is = new string("#INCLUSive");
  
  if (type == READ)
  {

    // set writing parameters off
    _isOpenWriting = false;

    // open file for reading 
    _ifs.open(fileName->c_str(), ios::in);
    if (_ifs)
    {
      // read first line
      getline(_ifs, _pLine, '\n');
      // cerr << "DEBUGLOG: " << _pLine << endl;
      if (_pLine.compare(0, 10, *is) == 0)
      {
        // this is a correct matrix file to start with
        _isOpenReading = true;
      }
      else
      {
        // something is wrong with this matrix file
        cerr << "--Error-- PMWIO(): Uncorrect format of matrix file: " <<
          *fileName << endl;
        cerrstr << "--Error--PMWIO(): Uncorrect format of matrix file: " <<
          *fileName << endl;
        cerrstr << "--Error--PMWIO(): First line reads: \n" << _pLine <<
          " (file should start with '#INCLUSive')" << endl;
        _pError = new string(cerrstr.str());
        cerrstr.flush();
        _isOpenReading = false;
        _ifs.close();
      }
    }
    else
    {
      cerr << "--Error-- PMWIO(): Unable to open file: " << *fileName <<
        " for reading." << endl;
      cerrstr << "--Error-- PMWIO(): Unable to open file: " << *fileName <<
        " for reading." << endl;
      _pError = new string(cerrstr.str());
      cerrstr.flush();
      _isOpenReading = false;
    }
  }
  else if (type == WRITE)
  {

    // set reading parameters off
    _isOpenReading = false;

    // open file for writing
    _ofs.open(fileName->c_str(), ios::trunc);
    if (_ofs)
    {
      _isOpenWriting = true;

      // write header line
      _ofs << "#INCLUSive Motif Model" << endl << "#" << endl;
      _ofs << "#The results in this file are NOT complete as long as there is no ";
      _ofs << "END-sign at the bottom of this file!" << endl;

    }
    else
    {
      cerr << "--Error-- PMWIO(): Unable to open file " << *fileName <<
        " for writing." << endl;
      cerrstr << "--Error-- PMWIO(): Unable to open file " << *fileName <<
        " for writing." << endl;
      _pError = new string(cerrstr.str());
      cerrstr.flush();
      _isOpenWriting = false;
    }
  }
  else
  {
    cerr << "--Error-- PMWIO(): Uncorrect type to open file " << *fileName <<
      " Use READ=0 || WRITE=1" << endl;
      cerrstr << "--Error-- PMWIO(): Uncorrect type to open file " << *fileName <<
      " Use READ=0 || WRITE=1" << endl;
      _pError = new string(cerrstr.str());
      cerrstr.flush();
    _isOpenWriting = false;
    _isOpenReading = false;
  }
  
  delete is;
  
}


/******************************************************************************
  Method:       delete
  Class:        PWMIO
  Arguments:    none
  
  Description:  destructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
PWMIO::~PWMIO()
{
  if (_pError != NULL) delete _pError; _pError = NULL;

  if (_isOpenReading && _ifs) {_ifs.close();}

  if (_isOpenWriting && _ofs)
	{
	_ofs <<"#END - Application ends." << endl;
    _ofs.close();
	}

	_pLine.clear();
	
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
  Method:       ReadMatrix
  Class:        PWMIO
  Arguments:    none
  
  Description:  get the next matrix from the matrix file
								a new PWM object is created 

  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  Description:  added reading of weight (MotifComparison BLiC)
  Author_Date:  MC_2009/07/21  
******************************************************************************/
PWM *
PWMIO::ReadMatrix()
{ 
  int W = 0, 
    nc = 0;
  float fA = 0.25,
    fC = 0.25,
    fG = 0.25,
    fT = 0.25;
  char * pcTag = NULL, 
    * pcValue = NULL;
  string * sId = NULL,
    * sCons = NULL;
  double score = 0;
  double weight(1); // 2009/07/19
  double *snf = 0;
  double **pM = NULL;
  bool bW = false,
    bCons = false,
    bID = false;
  PWM * pModel = NULL;
  ostringstream cerrstr;
  if (_pError != NULL) delete _pError; _pError = NULL;

  if (!_isOpenReading || !_ifs)
  {
    cerr << "--Error-- PWMIO::ReadMatrix(): Trying to read from file which is not open for reading\n";
    cerrstr <<  "--Error-- PWMIO::ReadMatrix(): Trying to read from file which is not open for reading\n";
    _pError = new string(cerrstr.str());
    cerrstr.flush();
    _isOpenReading = 0;
    return NULL;
  }

  if (_ifs.eof())
  {
    _isOpenReading = 0;
    return NULL;
  }

  
  // define strings
  pcTag = new char[64];
  pcValue = new char[512];

  // define default snf
  snf = new double[4];
  snf[0] = 0.0001;
  snf[1] = 0.0001;
  snf[2] = 0.0001;
  snf[3] = 0.0001;

  while ( _pLine.empty() || _pLine.length() <= 1 || _pLine[0] == '#')
  {
    nc = sscanf(_pLine.c_str(), "%s = %s", pcTag, pcValue);
    if (nc == 2)
    {
      if (strcmp(pcTag, "#ID") == 0)
      {
        sId = new string(pcValue);
        bID = true;
      }
      else if (strcmp(pcTag, "#Consensus") == 0)
      {
        sCons = new string(pcValue);
        bCons = true;
      }
      else if (strcmp(pcTag, "#W") == 0)
      {
        W = atoi(pcValue);
        // define matrix to store 
        if ( W > 0 )
        {
          bW = true;
          pM = new double*[W];
          for (int i = 0; i < W; i++)
            pM[i] = new double[4];
        }
        else
        {
          bW = false;
        }
      }
      else if (strcmp(pcTag, "#Weight") == 0) // 
      { 
        weight = (double) atof(pcValue);
      }
      else if (strcmp(pcTag, "#Score") == 0)
      {
        score = (double) atof(pcValue);
      }
    }

    if (!_ifs.eof())
    {
      getline(_ifs, _pLine, '\n');
    }
    else
    {
      _isOpenReading = 0;
      bW = false;
      break;
    }
  }
  
  // clean up temporary variables
  delete[] pcTag;
  delete[] pcValue;

  if ( bW )
  { 
    int i = 0;
    for (; i < W; i++)
    {
      // check the format of the input line
      int nc = sscanf(_pLine.c_str(), "%f%f%f%f", &fA, &fC, &fG, &fT);
      if (nc == 4)
      {
        pM[i][0] = (double)fA;
        pM[i][1] = (double)fC;
        pM[i][2] = (double)fG;
        pM[i][3] = (double)fT;
      }
      else
      {
        //cerr << "--Error-- PWMIO::ReadMatrix(): line is not as expected: " <<
          //_pLine << " in PWM with id=" << *sId << endl;
        cerrstr <<  "--Error-- PWMIO::ReadMatrix(): line (" <<
          _pLine << ") is not as expected in PWM";
        if (bID) cerrstr << " with ID=" << *sId;
        cerrstr << endl;
        _pError = new string(cerrstr.str());
        cerrstr.flush();
        break;
      }
      
      if (!_ifs.eof())
      {
        getline(_ifs, _pLine, '\n');
      }
      else
      {
        _isOpenReading = 0;
        for (int i = 0; i < W; i++)
          delete[] pM[i];
        delete[] pM;
        bW = false;
        break;
      }
    }
    if (i < W){ bW = false;} // reset so that no PWM is created
  }

  if ( bW )
  {
    // create PWM object
    pModel = new PWM(W, pM, snf);
    if ( !bID )
      sId = new string("unknown");
    pModel->SetID(sId);
    if ( bCons )
      pModel->SetConsensus(sCons);
    
    pModel->SetScore(score);
    pModel->SetWeight(weight); // 

    // delete local variables
    for (int i = 0; i < W; i++)
      delete[] pM[i];
    delete[] pM;
  }
  else
  {
    pModel = NULL; 
  }

  delete[] snf;

  if ( sId != NULL )
    delete sId;
  if ( sCons != NULL )
    delete sCons;
  return pModel;
}


/******************************************************************************
  Method:       WriteMatrix
  Class:        PWMIO
  Arguments:    PWM * pMatrix
  
  Description:  writes the PWM in the opened matrix file
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
  PWMIO::WriteMatrix(PWM * pMatrix)
{
  int W = 0,
    i = 0,
    j = 0;
  string *pId = NULL,
    * pConsensus = NULL;
  ostringstream cerrstr;
  if (_pError != NULL) delete _pError; _pError = NULL;

  if (_ofs && _isOpenWriting)
  {
    // write motif information
    W = pMatrix->Length();
    pId = pMatrix->GetID();
    pConsensus = pMatrix->GetConsensus();

    // cerr << "matrix: " << *pId << "    " << *pConsensus << endl;
    
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
      cerrstr <<  "--Error-- PWMIO::WriteMatrix(): File not open for writing" <<
      endl;
      _pError = new string(cerrstr.str());
      cerrstr.flush();
    return -1;
  }
}



/******************************************************************************
  Method:       IsOpen
  Class:        PWMIO
  Arguments:    none
  
  Description:  indicates if a file is open
  
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

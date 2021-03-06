#include "BackgroundIO.h"

#include <stdio.h>
#include <math.h>

// 2011/05/16 change ios::out into ios::trunc (first empty the outputfile)

/******************************************************************************
  Method:       BackgroundIO
  Class:        BackgroundIO
  Arguments:    string fileName
                int type (READ or WRITE)
  Description:  class constructor
                opens a file for reading or writing an INCLUSive background
                model. 
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
BackgroundIO::BackgroundIO(string fileName, int type)
{
  _pError = NULL;
  ostringstream cerrstr;
  if (_pError != NULL) delete _pError; _pError = NULL;
  if (type == READ)
  {

    _isOpenWriting = false;

    _ifs.open(fileName.c_str(), ios::in);

    if (_ifs)
    {
      _isOpenReading = true;

      // read first line
      getline(_ifs, _pLine, '\n');
      string is = "#INCLUSive";
      if (_pLine.compare(0, 10, is) == 0)
      {
        // this is a correct matrix file to start with
        _isOpenReading = true;
      }
      else
      { 
        // something is wrong with this matrix file
        cerr << "--Error-- BackgroundIO(): Uncorrect format of input file: "
          << fileName << endl;
        cerr << "--Error-- BackgroundIO(): First line reads: \n" << _pLine <<
          " (file should start with '#INCLUSive')" << endl;
        cerrstr << "--Error-- BackgroundIO(): Uncorrect format of input file: "
          << fileName << endl;
        cerrstr << "--Error-- BackgroundIO(): First line reads: \n" << _pLine <<
          " (file should start with '#INCLUSive')" << endl;
        _pError = new string(cerrstr.str());
        cerrstr.flush();
        _isOpenReading = false;
        _ifs.close();
      }
    }
    else
    {
      cerr << "--Error-- BackgroundIO::BackgroundIO() unable to open file: "
        << fileName << " for reading." << endl;
        cerrstr <<"--Error-- BackgroundIO::BackgroundIO() unable to open file: "
        << fileName << " for reading." << endl;
        _pError = new string(cerrstr.str());
        cerrstr.flush();
      _isOpenReading = false;
    }

  }
  else if (type == WRITE)
  {
    // set reading parameters off
    _isOpenReading = false;

    _ofs.open(fileName.c_str(), ios::trunc);
    if (_ofs)
    {
      _isOpenWriting = true;

      // write header line
      _ofs << "#INCLUSive Background Model v1.0" << endl << "#" << endl;
    _ofs << "#The results in this file are NOT complete as long as there is no ";
    _ofs << "END-sign at the bottom of this file!" << endl;

    }
    else
    {
      cerr << "--Error-- BackgroundIO(): Unable to open file for writing: " <<
        fileName << endl;
        cerrstr <<"--Error-- BackgroundIO(): Unable to open file for writing: " <<
        fileName << endl;
        _pError = new string(cerrstr.str());
        cerrstr.flush();
      _isOpenWriting = false;
    }

  }
  else
  {
    cerr << "--Error-- BackgroundIO(): Uncorrect type to open file "
      << fileName << ".  Use READ=0 -- WRITE=1\n" << endl;
        cerrstr << "--Error-- BackgroundIO(): Uncorrect type to open file "
      << fileName << ".  Use READ=0 -- WRITE=1\n" << endl;
        _pError = new string(cerrstr.str());
        cerrstr.flush();
    _isOpenWriting = false;
    _isOpenReading = false;
  }
}


/******************************************************************************
  Method:       ~BackgroundIO()
  Class:        BackgroundIO        
  Arguments:    none
  
  Description:  destructor
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
BackgroundIO::~BackgroundIO()
{
  if (_isOpenReading && _ifs)
    _ifs.close();

  if (_isOpenWriting && _ofs)
    _ofs.close();

  if (_pError != NULL) delete _pError; _pError = NULL;

  return;
}


/******************************************************************************
  Method:       Close
  Class:        BackgroundIO        
  Arguments:    none
  
  Description:  Close open file handle
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
void
BackgroundIO::Close()
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
  Method:       ReadBackgroundModel
  Class:        BackgroundIO        
  Arguments:    none
  
  Description:  Read background model from open file handle
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
BackgroundModel *
BackgroundIO::ReadBackgroundModel()
{

  int order = 0,
    i = 0,
    l = 0;
  string * organism = NULL;
  float fA, fC, fG, fT;
  char *pcTag, *pcValue;
  double * pSnf = NULL, *pOligo = NULL;
  double ** pTrans = NULL;
  bool bOrder = false,
    bGood = true;
  BackgroundModel *pBg = NULL;
  ostringstream cerrstr; 
  if (_pError != NULL) delete _pError; _pError = NULL;
  
	// variables to store tag and value in bg file
  pcTag = new char[16];
  pcValue = new char[256];

  if (!_isOpenReading || !_ifs)
  {
    cerr << "--Error-- BackgroundIO::ReadMatrix: Trying to read from file which is not open for reading.\n";
    cerrstr << "--Error-- BackgroundIO::ReadMatrix: Trying to read from file which is not open for reading.\n";
    _pError = new string(cerrstr.str());
    cerrstr.flush();
    // free some local variables
    delete[] pcValue;
    delete[] pcTag;

    return NULL;
  }

  while (bGood && !_ifs.eof())
  {
    // read next line
    getline(_ifs, _pLine, '\n');

    // tag lines
    if (!_pLine.empty() && sscanf(_pLine.c_str(), "%s = %s", pcTag, pcValue) == 2)
    {
      if (strcmp(pcTag, "#Order") == 0)
      {
        // order of background model
        order = atoi(pcValue);
        bOrder = true;
        
        // create internal variables to store snf, transition matrix and oligo frequencies
        l = (int) pow(4.0, order);
        if ( order == 0 )
          l = 4;  // 
        pSnf = new double[4]; // single nucleotide frequency
        pOligo = new double[l]; // oligo frequencies
        pTrans = new double*[l];
        for (i=0; i<l; i++)
          pTrans[i] = new double[4];

      }
      else if (strcmp(pcTag, "#Organism") == 0)
      {
        // organism name
        organism = new string(pcValue);
      }
    }
    
    // data lines
    if (bOrder && !_pLine.empty() && strncmp(_pLine.c_str(), "#snf", 4) == 0)
    {
      // this the snf part -> next line should contain the values
      getline(_ifs, _pLine, '\n');

      // check the format
      if (sscanf(_pLine.c_str(), "%f%f%f%f", &fA, &fC, &fG, &fT) == 4)
      {
        pSnf[0] = (double) fA;
        pSnf[1] = (double) fC;
        pSnf[2] = (double) fG;
        pSnf[3] = (double) fT;
      }
      else
      {
        cerr <<
          "--Error-- BackgroundIO::ReadBackgroundModel(): Uncorrect format of the input file while trying to extract SNF line:\n  "
          << _pLine << endl;
    cerrstr << "--Error-- BackgroundIO::ReadBackgroundModel(): Uncorrect format of the input file while trying to extract SNF line:\n  "
          << _pLine << endl;
    _pError = new string(cerrstr.str());
    cerrstr.flush();
        bGood = false;
        break;
      }
    }
    else if (bOrder && strncmp(_pLine.c_str(), "#oligo", 6) == 0)
    {
      // this is the set of oligo frequencies
      // next 4^Order lines make up the oligo frequencies
      for (i=0; i<l; i++)
      {
        // read next line 
        getline(_ifs, _pLine, '\n');

        if (sscanf(_pLine.c_str(), "%f", &fA) == 1)
        {
          // correct format
          pOligo[i] = (double) fA;
        }
        else
        {
          // wrong format -> return empty pointer
          cerr << "--Error-- BackgroundIO::ReadBackgroundModel(): Uncorrect format of the input file while trying to extract oligo frequency line:\n" << _pLine << endl;
    cerrstr << "--Error-- BackgroundIO::ReadBackgroundModel(): Uncorrect format of the input file while trying to extract oligo frequency line:\n" << _pLine << endl;
    _pError = new string(cerrstr.str());
    cerrstr.flush();
          bGood = false;
          break;
        }
      }
    }
    else if (bOrder && !_pLine.empty() && strncmp(_pLine.c_str(), "#transition", 11) == 0)
    {
      // this is the transition matrix
      
      for (i=0; i<l; i++)
      {
        // get next line
        getline(_ifs, _pLine, '\n');

        if (sscanf(_pLine.c_str(), "%f %f %f %f", &fA, &fC, &fG, &fT) == 4)
        {
          // correct format
          pTrans[i][0] = (double) fA;
          pTrans[i][1] = (double) fC;
          pTrans[i][2] = (double) fG;
          pTrans[i][3] = (double) fT;

        }
        else
        {
          cerr << "Error\n BackgroundIO::ReadBackgroundModel: Uncorrect format of the input file while trying to extract transition matrix line:\n   "
            << _pLine << endl;
    cerrstr <<  "Error\n BackgroundIO::ReadBackgroundModel: Uncorrect format of the input file while trying to extract transition matrix line:\n   "<< _pLine << endl;
    _pError = new string(cerrstr.str());
    cerrstr.flush();
          bGood = false;
          break;
        }
      }
    }
    else
    {
      continue;
    }
  }

  // create the background model
  if ( bGood )
  {
    pBg = new BackgroundModel(order, pTrans, pOligo, pSnf);
    if (organism != NULL && !organism->empty())
      pBg->SetOrganism(organism);
  } 
  else
  {
    pBg = NULL;
  }
  
  // free some local variables
  if (pcValue != NULL)
    delete[]pcValue;
  if (pcTag != NULL)
    delete[]pcTag;
  if (pSnf != NULL)
    delete[]pSnf;
  if (pOligo != NULL)
    delete[]pOligo;
  if (pTrans != NULL)
  {
    for (i=0; i<l; i++)
      delete[] pTrans[i];
    delete[] pTrans;
  }
  
  return pBg;

}


/******************************************************************************
  Method:       Isopen()
  Class:        BackgroundIO        
  Arguments:    none
  
  Description:  Method to check if the file is open for reading or writing
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
bool
BackgroundIO::IsOpen()
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


/******************************************************************************
  Method:       WriteBackgroundModel
  Class:        BackgroundIO        
  Arguments:    BackgroundModel *bgM
  
  Description:  Method to write background model in INCLUSive format
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
int
BackgroundIO::WriteBackgroundModel(BackgroundModel * bgM)
{
  ostringstream cerrstr; 
  if (_pError != NULL) delete _pError; _pError = NULL;
  int i, j, order, W;

  if (_ofs && _isOpenWriting)
  {
    // print header
    order = bgM->GetOrder();
    _ofs << "#Order = " << order << endl;
    if (bgM->GetOrganism() != NULL)
    {
      _ofs << "#Organism = " << *(bgM->GetOrganism()) << endl;
    }
    else
    {
      _ofs << "#Organism = unknown" << endl;
    }
    if (bgM->GetFileName() != NULL)
    {
      _ofs << "#Sequences = " << *(bgM->GetFileName()) << endl;
    }
    else
    {
      _ofs << "#Sequences = not_defined" << endl;
    }
    //_ofs << "#Path = " << endl;
    _ofs << "#" << endl << endl;

    // print snf
    _ofs << "#snf" << endl;
    for (i = 0; i < 4; i++)
    {
      _ofs << bgM->GetSnfValueAt(i) << "\t";
    }
    _ofs << endl << endl;

    // print oligo frequencies
    _ofs << "#oligo frequency" << endl;
    W = (int) pow(4.0, order);
		if ( order == 0 )
			W = 4;
    for (i = 0; i < W; i++)
    {
      _ofs << bgM->GetOligoFrequencyValueAt(i) << endl;
    }
    _ofs << endl;

    // print transition matrix
    _ofs << "#transition matrix" << endl;
    for (i = 0; i < W; i++)
    {
      for (j = 0; j < 4; j++)
      {
        _ofs << bgM->GetTransitionMatrixValueAt(i, j) << "\t";
      }
      _ofs << endl;
    }
    _ofs << endl;
    _ofs << "#END - BackgroundModel ends." << endl;

  }
  else
  {
    cerr <<
      "--Error-- BackgroundIO::WriteBackgroundModel(): unable to write to file."
      << endl;
    cerrstr <<  "--Error-- BackgroundIO::WriteBackgroundModel(): unable to write to file."
      << endl;
    _pError = new string(cerrstr.str());
    cerrstr.flush();
    return 0;
  }

  return 1;
}

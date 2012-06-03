#include "GFFWriter.h"
#include <iostream>


// 2011/05/16 change ios::out into ios::trunc (first empty the outputfile)

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
GFFWriter::GFFWriter()
{
  // write to STDOUT
  _isOpenWriting = true;
  _isSTDOUT = true;
  cout << "#INCLUSive GFF File" << endl;
}

/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
GFFWriter::GFFWriter(string * pFileName)
{
  // open file
  _ofs.open(pFileName->c_str(), ios::trunc);
  if (_ofs)
  {
    _isOpenWriting = true;

    // write header line
    _ofs << "#INCLUSive GFF File" << endl;
    _ofs << "#" << endl;
    _ofs << "#The results in this file are NOT complete as long as there is no ";
    _ofs << "END-sign at the bottom of this file!" << endl;

  }
  else
  {
    cerr << "--Error-- GFFWriter(): Unable to open file " << *pFileName <<
      " for writing." << endl;
    _isOpenWriting = false;
  }

  _isSTDOUT = false;

}


/******************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
******************************************************************************/
GFFWriter::~GFFWriter()
{
  if (_isOpenWriting && !_isSTDOUT)
	{
    _ofs << "#END - Application ends." << endl;
    _ofs.close();
	}

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
  GFFWriter::WriteInstance(Instance * pInst, string * source)
{
  int
    length = pInst->ParentSequence()->Length();

  if (_isSTDOUT)
  {
    cout << *(pInst->ParentSequence()->GetID()) << "\t";
    cout << *source << "\t";
    cout << "misc_feature\t";

    if (pInst->Strand() == plus_strand)
    {
      cout << pInst->Start() + 1 << "\t";
      cout << pInst->Start() + pInst->Length() << "\t";
      cout << pInst->Score() << "\t";
      cout << "+\t";
    }
    else
    {
      cout << length - pInst->Start() - pInst->Length() + 1 << "\t";
      cout << length - pInst->Start() << "\t";
      cout << pInst->Score() << "\t";
      cout << "-\t";
    }
    cout << ".\t";
    cout << "id \"" << *(pInst->ID()) << "\"; ";
    cout << "site \"" << *(pInst->PrintSite()) << "\"; " << endl;
  }
  else
  {
    _ofs << *(pInst->ParentSequence()->GetID()) << "\t";
    _ofs << *source << "\t";
    _ofs << "misc_feature\t";
    if (pInst->Strand() == plus_strand)
    {
      _ofs << pInst->Start() + 1 << "\t";
      _ofs << pInst->Start() + pInst->Length() << "\t";
      _ofs << pInst->Score() << "\t";
      _ofs << "+\t";
    }
    else
    {
      _ofs << length - pInst->Start() - pInst->Length() + 1 << "\t";
      _ofs << length - pInst->Start() << "\t";
      _ofs << pInst->Score() << "\t";
      _ofs << "-\t";
    }
    _ofs << ".\t";
    _ofs << "id \"" << *(pInst->ID()) << "\"; ";
    _ofs << "site \"" << *(pInst->PrintSite()) << "\"; " << endl;
  }

  return;
}


void
GFFWriter::AddComment(string * pComStr)
{
  if (pComStr == NULL)
    return;

  if (_isSTDOUT)
  {
    cout << *pComStr << endl;
  }
  else
  {
    _ofs << *pComStr << endl;
  }
}
/******************************************************************************
  Description:  GFFWriter::AddComment(const string &input)
  Author_Date:  MC_2009/01/15
******************************************************************************/
void GFFWriter::AddComment(const string &input)
{
  if (_isSTDOUT) cout << input << endl;
  else _ofs << input << endl;
}


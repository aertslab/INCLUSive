#include "GFFWriter.h"
#include <iostream>


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
  _ofs.open(pFileName->c_str(), ios::out);
  if (_ofs)
  {
    _isOpenWriting = true;

    // write header line
    _ofs << "#INCLUSive GFF File" << endl;

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
    _ofs.close();

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

#include "GFFReader.h"

GFFReader::GFFReader(string * gffFile )
{
  // open file for reading
  _ifs.open(gffFile->c_str(), ios::in);
  
  if (!_ifs)
  {
    // file is not open for reading
    _bIsOpenReading = false;
    cerr << "GFFReader::GFFReader(): Unable to open " << *gffFile << " for reading." << endl;
  }
  else
  {
    _bIsOpenReading = true;
  }
}


GFFReader::~GFFReader()
{
  if ( _pLine != NULL )
    delete _pLine;
  _pLine = NULL;
  
  if ( _ifs.good())
    _ifs.close();

  return;
}


Instance * 
GFFReader::NextInstance()
{
  Instance * pFeat = NULL;
  
  do // read input as long as a # is encountered => comment lines
  {
    getline(_ifs, _pLine, '\n');   
  }
  while ( ( _pLine[0] != '#' || _ifs.empty() ) && !_ifs.eof());

  // parse the GFF line
  string tokens(" \t");
  string::size_type pos = 0;

  if ((pos = _pLine->find_first_of(tokens)) != string::npos)
  {
    
  }
  else
  {
    
  }
}

#include "GFFDrawer.h"
#include <time.h>

const int GFFDrawer::_xOffSet = 20;
const int GFFDrawer::_yOffSet = 20;
const int GFFDrawer::_seqOffSet = 40;

GFFDrawer::GFFDrawer(string *pFileName, list<SequenceObject*>* pSeqSet)
{
  // open file for writing
  _ofs.open(pFileName->c_str(), ios::out);
  if (_ofs)
  {
    // file is open for writing
    _isOpenWriting = true;
    _pFile = pFileName;
    
    // define dimensions of the plot
    _nbrSequences = pSeqSet->size();
    _maxLength = 0;
    
    list<SequenceObject*>::iterator si = pSeqSet->begin();
    _nbrSequences = 0;
    while ( si != pSeqSet->end())
    {
      // store pointers to sequence object in a map
      _seqMap.insert(map<int, SequenceObject*>::value_type(_nbrSequences, (*si)));
      _nbrSequences++;
      
      // check length
      int l = (*si)->Length();
      _maxLength = ( _maxLength < l ) ? l : _maxLength;
      si++;
    }
  
    // define picture dimensions
    _xDim = _xOffSet + _maxLength + _xOffSet;
    _yDim = _yOffSet + ( _nbrSequences * _seqOffSet) + _yOffSet;

    // initialize figure with EPS headers
    _InitiateFigure();

    // draw basic figure outline
    _DrawOutline();
    
  }
  else
  {
    cerr << "--Error-- GFFDrawer(): Unable to open file " << *pFileName <<
      " for writing." << endl;
    _isOpenWriting = false;
  }

}


GFFDrawer::~GFFDrawer()
{
  // 
}

void 
GFFDrawer::AddInstance(Instance* PInst)
{

  return;
}


void 
GFFDrawer::Close()
{
  _CloseFigure();
  return;
}


// local functions
void 
  GFFDrawer::_InitiateFigure()
{
  // char pDate[128];

  // get the date 
  time_t rtime;
  time(&rtime);
  // strftime(pDate, 128, "%a %b %d %H:%M:%S %Y", &rtime);
  
  _ofs << "%!PS-Adobe-3.0 EPSF-3.0"  << endl
      << "%%BoundingBox: 0 0 " << _xDim << " " << _yDim << endl
      << "%%Title: INCLUSive - Motif Alginment Plot"  << endl
      << "%%Creator: INCLUSive"  << endl
      << "%%CreationDate: " << ctime(&rtime) << endl
      << "%%DocumentFonts: Helvetica Helvetica-Bold"  << endl
      << "%%DocumentNeededFonts: Helvetica Helvetica-Bold"  << endl
      << "%%Pages: 1"  << endl
      << "%%EndComments" << endl;
 
  // define eps fonts
  _ofs << "/LABELfont"  << endl
    << "{"  << endl
    << "  /Helvetica-Bold findfont 10 scalefont setfont"  << endl
    << "}def"  << endl;
 
 return;
}

void
  GFFDrawer::_DrawOutline()
{
  // draw ruler
  _ofs << "newpath " <<  endl
    << _xOffSet << " " << (_yOffSet + (_seqOffSet * _nbrSequences)) << " moveto" << endl
    << _maxLength << " " << 0 << " rlineto" << endl
    << "2 setlinewidth" << endl
    << "stroke" <<  endl;

  // draw start axis 
  _ofs << "newpath" <<  endl
    << (_xOffSet + _maxLength) << " " << _yOffSet << " moveto" << endl
    << "0" << " " << (_seqOffSet * _nbrSequences) << " rlineto" <<endl
    << "2 setlinewidth" << endl
    << "stroke" <<  endl;

  // draw thicks at every 50bp in the ruler 
  
  
  // draw sequences
  
  return;
}

void
  GFFDrawer::_CloseFigure()
{
  _ofs << "end" << endl
    << "showpage" << endl
    << "%%EOF" << endl;
  return;
}

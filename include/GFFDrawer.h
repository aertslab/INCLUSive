#ifndef gffdrawer_include_declared
#define gffdrawer_include_declared

#include "inclusive.h"
#include "utilities.h" 
#include "Instance.h"
#include "SequenceObject.h"

struct color 
{
  uint red;
  uint green;
  uint blue;
};

typedef pair<string*, color*> mColor;

class GFFDrawer 
{
  private:
    ofstream _ofs;
    uint _nbrSequences;
    int _maxLength;
    string* _pFile;
    bool _isOpenWriting;
  
    list<mColor *> _colorList;
    map<int, SequenceObject*> _seqMap;
  
  // dimension parameters
  int _xDim;
  int _yDim;
  static const int _xOffSet;
  static const int _yOffSet;
  static const int _seqOffSet;

  // drawing functions
  void _InitiateFigure();
  void _CloseFigure();
  void _AddSequenceLine();
  void _AddInstance();
  void _DrawOutline();

  public:
    GFFDrawer(string *fileName, list<SequenceObject*> *pSeqSet);
    ~GFFDrawer();
  
    void AddInstance(Instance* pInst);
  
    void Close();
  
};

#endif

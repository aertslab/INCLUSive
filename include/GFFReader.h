#ifndef gffreader_include_declared
#define gffreader_include_declared

#include "inclusive.h"
#include "Instance.h"

class GFFReader 
{
  private:
    ifstream _ifs;
    bool _bIsOpenreading;
    string * _pLine;

  public:
    GFFReader(string* gffFile);
    ~GFFReader();
    void Close();
    
    Instance * NextInstance();    
    bool IsOpen(){ return _bIsOpenreading;};
    
};

#endif

#ifndef instance_include_declared
#define instance_include_declared

#include "inclusive.h"
#include "SequenceObject.h"

class SequenceObject;

class Instance{
 private:
  int _start;
  int _length;
  strand_modes _strand;
  double _score;
  string *_pId;
  Sequence *_pSite;
  string *_pSiteString;
  SequenceObject * _pParentSequence;
  string *_pSeqName; 	
 
 
 public:
	 // constructors
  Instance(Instance &copyInstance);
  Instance(SequenceObject *pSeqObj, strand_modes s, int start, int length);
  // destructor
  ~Instance();

  // inspectors
  int Start(){ return _start; };
  int Length(){ return _length; };
  strand_modes Strand(){ return _strand;};
  double Score(){ return _score; };
  string *ID(){ return _pId; };
  Sequence * Site(){ return _pSite; };
  string * PrintSite(){ return _pSiteString; };
  SequenceObject * ParentSequence(){return _pParentSequence;};
  string * PrintSeqName(){return _pSeqName;};

  bool CheckLeftExtension(int pos);
  bool CheckRightExtension(int pos);
  
  // adaptors
  void SetID(string *pId){ _pId = pId; return;};
  Instance * SetScore(double score){ _score = score; return this;};
  Instance * SetSeqName(string *pName){_pSeqName = pName; return this;};
  Instance * ExtendLeft(int pos);
  Instance * ExtendRight(int pos);
    
};



#endif

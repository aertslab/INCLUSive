#ifndef gffwriter_include_defined
#define gffwriter_include_defined

#include "inclusive.h"
#include "Instance.h"
#include "SequenceObject.h"

class GFFWriter
{
  private:
	bool _isOpenWriting;
	ofstream _ofs;
  bool _isSTDOUT;
  
  public:
	// constructor
	GFFWriter(string *pFileName);  // writes to file
  GFFWriter();	                 // writes to STDOUT
	
	// destructor
	~GFFWriter();
	  
	//   
	void WriteInstance(Instance *pInst, string *source);
  void AddComment(string *pComStr);
  void AddComment(const string &input);//
};

#endif

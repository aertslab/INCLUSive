#ifndef maskvector_include_declared
#define maskvector_include_declared

#include <vector>
#include "inclusive.h"

class MaskVector {
  private:
	vector<int> * _mask;
  
  public:  
	// constructor  
	MaskVector(int l);
  
	// destructor
	~MaskVector();
  
  
	void ResetMask(); // make all one
	void ClearMask(); // make all zero
	void UpdateMask(int start, int w, int value); // set value from start to start+w-1
	void SetValueAt(int i, int value);
  
	// inspectors
	int GetValueAt(int i);
	vector<int> * GetMask(){ return _mask;};
	MaskVector * ReverseMask();
  
  vector<int>::iterator begin(){return _mask->begin();};
  vector<int>::iterator end(){return _mask->end();};
  
};

#endif

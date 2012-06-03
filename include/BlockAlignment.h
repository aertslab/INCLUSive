#ifndef blockalignment_declared
#define blockalignment_declared

#include "PWM.h"

class BlockAlignment
{
  private:
    double ** _pAlignmentMatrix;
    int _length1;
    int _length2;
    int _maxXindex; 
    int _maxYindex; 
    int _pathXStart;
    int _pathYStart;  
    double _maxValue;
  
  public:
    // constructor
    BlockAlignment(int w1, int w2);
  
    // destructor 
    ~BlockAlignment();
    
    // adaptors
    void UpdateMatrix(PWM* pwm1, PWM* pwm2, double threshold, double gapscore);
    void ResetMatrix();
  
    // inspector
    int GetMaxXIndex(){return _maxXindex;};
    int GetMaxYIndex(){return _maxYindex;};
    double GetMaxValue(){return _maxValue;};
    double GetValueAt(int i, int j);
    int GetPathLength(){return _maxXindex - _pathXStart + 1;};
    void GetPathStart(int* x, int* y);
    
};

#endif

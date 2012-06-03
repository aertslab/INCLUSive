#include "BlockAlignment.h"
#include <math.h>

BlockAlignment::BlockAlignment(int w1, int w2)
{

  _length1 = w1;
  _length2 = w2;

  // create internal alignment matrices
  _pAlignmentMatrix = new double*[w1];
  for ( int i=0; i<_length1; i++)
  {
    _pAlignmentMatrix[i] = new double[_length2];
    for ( int j=0; j<_length2; j++)
      _pAlignmentMatrix[i][j] = 0;
  }
  // initialize internal variables
  _maxXindex = 0;
  _maxYindex = 0;
  _maxValue = 0;
  _pathXStart = 0;
  _pathYStart = 0;

}


BlockAlignment::~BlockAlignment()
{
  for ( int i=0; i<_length1; i++)
    delete[] _pAlignmentMatrix[i];
  delete[] _pAlignmentMatrix;
  
  _length1 = 0;
  _length2 = 0;
}


void BlockAlignment::ResetMatrix()
{
  for ( int i=0; i<_length1; i++)
    for ( int j=0; j<_length2; j++)
      _pAlignmentMatrix[i][j] = 0;
  return;
}
  

void 
BlockAlignment::UpdateMatrix(PWM* pwm1, PWM* pwm2, double threshold, double gapScore)
{
  double score = 0, L1 = 0, L2 = 0;
  
  _maxValue = 0;
  _maxXindex = 0;
  _maxYindex = 0;
  
  for (int i=0; i<_length1; i++)
  {
    for (int j=0; j<_length2; j++ )
    {
      score = threshold;
      for (int k=0; k<4; k++)
        score -= ((pwm1->GetValueAt(i,k) -  pwm2->GetValueAt(j,k)) * (log(pwm1->GetValueAt(i,k)) - log(pwm2->GetValueAt(j,k))));
      
      if ( i && j )
      {
        L1 = _pAlignmentMatrix[i-1][j-1] + score;
        L2 = _pAlignmentMatrix[i-1][j-1] - gapScore;
      }
      else
      {             
        L1 = score;
        L2 = -gapScore;
      }
      
      if ( L1 < 0 && L2 < 0 )
      {
        _pAlignmentMatrix[i][j] = 0;
      }
      else if ( L1 > L2 )
      {
        _pAlignmentMatrix[i][j] = L1;
      }
      else
      {
        _pAlignmentMatrix[i][j] = L2;
      }

      if ( _pAlignmentMatrix[i][j] > _maxValue )
      { 
        _maxValue = _pAlignmentMatrix[i][j];
        _maxXindex = i;
        _maxYindex = j;
      }
    }
  }  
  
  // compute path start and length
  _pathXStart = _maxXindex;
  _pathYStart = _maxYindex;
  while ( _pathXStart && _pathYStart && _pAlignmentMatrix[_pathXStart-1][_pathYStart-1] > 0 )
  {
    _pathXStart--;
    _pathYStart--;
  }

  return;
}



void 
BlockAlignment::GetPathStart(int* x, int* y)
{
  *x = _pathXStart;
  *y = _pathYStart;
  return;
}

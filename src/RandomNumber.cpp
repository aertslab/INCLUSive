#include "RandomNumber.h"
#include <iostream>
#include <time.h>
#include <cmath>
// release version 3.1.2 : 1 bug fixed in GetUniform

// initialization of the internal state of the random number generator
long RandomNumber::_lrand1 = 1000;
long RandomNumber::_lrand2 = 2000;
long RandomNumber::_lrand3 = 3000;

/*************************************************************************
  Method:       InitRandomNumber
  Class:        RandomNumber
  Arguments:    none
  
  Description:  initialization of the random number generator
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*************************************************************************/
void
RandomNumber::InitRandomNumber()
{
  _lrand1 = 1000 * (long) time(NULL);
  _lrand2 = 2000 * (long) time(NULL);
  _lrand3 = 3000 * (long) time(NULL);
} 


/*************************************************************************
  Method:       _lcg
  Class:        RandomNumber
  Arguments:    int a, int c, long M, long R
  
  Description:  internal function of the random number generator
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*************************************************************************/
long
RandomNumber::_lcg(int a, int c, long M, long R)
{
  return ((a * R + c) % M);
}

/*************************************************************************
  Method:       _mod
  Class:        RandomNumber
  Arguments:    float a, float b
  
  Description:  modulo(a,b)
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*************************************************************************/
float
RandomNumber::_mod(float a, float b)
{
  return a - b * floor(a / b);
}


/*************************************************************************
  Method:       GetUniform
  Class:        RandomNumber
  Arguments:    none
  
  Description:  get a uniform distributed random number
  
  Date:         2005/12/06
  Author:       Marleen Claeys
  
*************************************************************************/
float
RandomNumber::GetUniform()
{
  // update internal state
  _lrand1 = _lcg(171, 0, RAND1, _lrand1);  
  _lrand2 = _lcg(172, 0, RAND2, _lrand2);  
  _lrand3 = _lcg(170, 0, RAND3, _lrand3);
    
    // version 3.1.2 : _lrand2 ipv _lrand3 in teller van RAND2
    return _mod((float) (_lrand1) / RAND1 + (float) (_lrand2) / RAND2 +
               (float) (_lrand3) / RAND3, 1);
} 


/*************************************************************************
  Method:       GetExponential
  Class:        RandomNumber
  Arguments:    float lambda
  
  Description:  
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*************************************************************************/
float
RandomNumber::GetExponential(float lambda)
{
  return -log(GetUniform()) / lambda;
}

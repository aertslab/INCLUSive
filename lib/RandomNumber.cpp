#include "RandomNumber.h"
#include <iostream>
#include <time.h>
#include <cmath>
// #include <libio.h>

// variable initialization
long RandomNumber::_lrand1 = 1000;
long RandomNumber::_lrand2 = 2000;
long RandomNumber::_lrand3 = 3000;

/*************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*************************************************************************/void
RandomNumber::InitRandomNumber()
{
  _lrand1 = 1000 * (long) time(NULL);
  _lrand2 = 2000 * (long) time(NULL);
  _lrand3 = 3000 * (long) time(NULL);
} 


/*************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*************************************************************************/long
RandomNumber::_lcg(int a, int c, long M, long R)
{
  return ((a * R + c) % M);
}

  /*************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*************************************************************************/
float
RandomNumber::_mod(float a, float b)
{
  return a - b * floor(a / b);
}


/*************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*************************************************************************/
float
RandomNumber::GetUniform()
{
  // update values
  _lrand1 = _lcg(171, 0, RAND1, _lrand1);  _lrand2 = _lcg(172, 0, RAND2, _lrand2);  _lrand3 = _lcg(170, 0, RAND3, _lrand3);
    return _mod((float) (_lrand1) / RAND1 + (float) (_lrand3) / RAND2 +
               (float) (_lrand3) / RAND3, 1);
} 


/*************************************************************************
  Method: 
  Class:        
  Arguments: 
  
  Description:
  
  Date:         2003/06/26
  Author:       Gert Thijs <gert.thijs@esat.kuleuven.ac.be>
  
*************************************************************************/
float
RandomNumber::GetExponential(float lambda)
{
  return -log(GetUniform()) / lambda;
}

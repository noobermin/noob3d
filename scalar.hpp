#ifndef _SCALAR_H_
#define _SCALAR_H_

#define NOOB3D_SMALL 1.0e-10

#include <cstdlib>
#include <cmath>

namespace noob3d
{
  typedef double scalar;

  inline scalar
  dist(scalar a, scalar b)
  {return abs(a-b);}

  inline bool
  eq(scalar l, scalar r)
  { return dist(l,r) <= NOOB3D_SMALL; }
  
  inline scalar
  inverse(scalar in)
  { return 1.0/in; }
  
  inline scalar
  squareRoot(scalar in)
  {return sqrt(in);}

  inline scalar
  sq(scalar in)
  {return in*in;}
}

#endif //_SCALAR_H_

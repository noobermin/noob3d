#ifndef _NOOB3D_VECTOR_H
#define _NOOB3D_VECTOR_H

#include <vector>
#include "scalar.hpp"
#include "scalar.hpp"

namespace noob3d
{

  ////////////////////
  // vector classes //
  ////////////////////
 
  class vector2d
  {
  public:
    scalar x, y;
    vector2d(scalar xin=0.0, scalar yin=0.0) 
      : x(xin), y(yin) {};

    vector2d(const vector2d& v)
      : x(v.x), y(v.y) {};
  
   /*public operators*/

    //accessors
    scalar
    operator[] (int index) const
    { return index ? y : x;}
    scalar&
    operator[] (int index)
    { return index ? y : x;}
    
    //normal arithmetic
    vector2d
    operator+ (const vector2d& rval) const
    {
      vector2d ret(*this);
      ret.x+=rval.x;
      ret.y+=rval.y;
      return ret;
    }
    vector2d
    operator- (const vector2d& rval) const
    {
      vector2d ret(*this);
      ret.x-=rval.x;
      ret.y-=rval.y;
      return ret;
    }

    vector2d
    operator* (scalar rval) const
    {
      vector2d ret(*this);
      ret.x*=rval;
      ret.y*=rval;
      return ret;
    }
    // vector
    //operator* (const matrix &rval) const;

    vector2d
    operator/ (scalar rval) const
    {
      vector2d ret(*this);
      ret.x/=rval;
      ret.y/=rval;
      return ret;
    }

    //More efficient self assignment
    vector2d&
    operator= (const vector2d& rval) //note, this is arithmetic assign, although there really is no diff...
    {
      x=rval.x; y=rval.y;
      return *this;
    }
    vector2d&
    operator+= (const vector2d& rval)
    {
      x+=rval.x; y+=rval.y;
      return *this;
    }

    vector2d&
    operator-= (const vector2d& rval)
    {
      x-=rval.x; y-=rval.y;
      return *this;
    }

    vector2d&
    operator*= (scalar rval)
    {
      x*=rval; y*=rval;
      return *this;
    }
    vector2d&
    operator/= (scalar rval)
    {
      x/=rval; y/=rval;
      return *this;
    }
       

    //comparators

    scalar
    squareLength() const
    { return x*x+y*y;}

    scalar
    length() const
    { return squareRoot(squareLength()); }
    
    vector2d&
    normalize()
    { return (*this)/=length(); }

    vector2d
    normalized() const
    { vector2d a(*this); return a.normalize();}

    bool
    operator== (vector2d &rval) const
    { return eq(((*this)-rval).length(),0.0);}
    bool
    operator!= (vector2d &rval) const
    { return !((*this)==rval); }
  };

  class vector3d
  {
  public:
    scalar x,y,z;
    vector3d(scalar xin=0.0, scalar yin=0.0, scalar zin=0.0) 
      : x(xin), y(yin), z(zin) {};

    vector3d(const vector3d& v)
      : x(v.x), y(v.y), z(v.z) {};
   /*public operators*/

    //accessors
    scalar
    operator[] (int index) const
    { return index ? (index == 2 ? z : y) : x;}
    scalar&
    operator[] (int index)
    { return index ? (index == 2 ? z : y) : x;}
    
    //normal arithmetic
    vector3d
    operator+ (const vector3d& rval) const
    {
      vector3d ret(*this);
      ret.x+=rval.x;
      ret.y+=rval.y;
      ret.z+=rval.z;
      return ret;
    }
    vector3d
    operator- (const vector3d& rval) const
    {
      vector3d ret(*this);
      ret.x-=rval.x;
      ret.y-=rval.y;
      ret.z-=rval.z;
      return ret;
    }

    vector3d
    operator* (scalar rval) const
    {
      vector3d ret(*this);
      ret.x*=rval;
      ret.y*=rval;
      ret.z*=rval;
      return ret;
    }

    vector3d
    operator/ (scalar rval) const
    {
      vector3d ret(*this);
      ret.x/=rval;
      ret.y/=rval;
      ret.z/=rval;
      return ret;
    }

    //More efficient self assignment
    vector3d&
    operator= (const vector3d& rval) //note, this is arithmetic assign, although there really is no diff...
    {
      x=rval.x; y=rval.y; z=rval.z;
      return *this;
    }
    vector3d&
    operator+= (const vector3d& rval)
    {
      x+=rval.x; y+=rval.y; z+=rval.z;
      return *this;
    }

    vector3d&
    operator-= (const vector3d &rval)
    {
      x-=rval.x; y-=rval.y; z-=rval.z;
      return *this;
    }

    vector3d&
    operator*= (scalar rval)
    {
      x*=rval; y*=rval; z*=rval;
      return *this;
    }
    vector3d&
    operator/= (scalar rval)
    {
      x/=rval; y/=rval; z/=rval;
      return *this;
    }

    scalar
    squareLength() const
    { return x*x+y*y+z*z;}

    scalar
    length() const
    { return squareRoot(squareLength()); }
    
    vector3d&
    normalize()
    { 
      return (*this)/=length();
    }

    vector3d
    normalized() const
    { vector3d a(*this); return a.normalize();}
    
    //comparators

    bool
    operator!= (vector3d& rval) const
    { return !((*this)==rval); }

    bool
    operator== (vector3d& rval) const
    { return eq(((*this)-rval).length(),0.0); }
  };


  inline scalar
  dot(const vector2d& lval, const vector2d& rval)
  {  
    return lval.x*rval.x + lval.y* rval.y;
  }

  inline scalar
  dot(const vector3d& lval, const vector3d& rval)
  {  
    return lval.x*rval.x + lval.y*rval.y + lval.z*rval.z;
  }

  inline vector3d
  cross(const vector3d& lval, const vector3d& rval)
  {
    vector3d ret;
    ret.x = lval.y * rval.z - lval.z * rval.y;
    ret.y = lval.z * rval.x - lval.x * rval.z;
    ret.z = lval.x * rval.y - lval.y * rval.x;
    return ret;
  }
  //returns only the z component (spin)
  inline scalar
  cross(const vector2d& lval, const vector2d& rval)
  {
    return lval[0] * rval[1] - lval[1] * rval[0];
  }

  inline vector2d
  normal( const vector2d& in)
  {
    return vector2d(in[1], -in[0]);
  }
  
  inline vector2d
  orthonormal(const vector2d& in)
  {
    return normal(in).normalize();
  }

  inline vector2d
  trunc3to2(const vector3d& in)
  {
    return vector2d(in[0], in[1]);
  }

  inline vector2d
  midpoint(const vector2d& r1, const vector2d& r2)
  {
    return (r1 + r2) / 2.0f;
  }
  inline vector3d
  midpoint(const vector3d& r1, const vector3d& r2)
  {
    return (r1 + r2) / 2.0f;
  }

  inline scalar
  dist(const vector2d& a, const vector2d& b)
  {return (a-b).length();}

  inline scalar
  dist(const vector3d& a, const vector3d& b)
  {return (a-b).length();}
}

#endif /* NOOB3D_VECTOR_H*/

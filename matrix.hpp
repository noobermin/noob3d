#ifndef _NOOB3D_MATRIX_H_
#define _NOOB3D_MATRIX_H_

#include <array>
#include <tuple>
#include <algorithm>
#include "scalar.hpp"
#include "vector.hpp"
#include "types.hpp"
#include <iostream>
#include <iomanip>

namespace noob3d
{
  class matrix2d
  {
    std::array<std::array<scalar,2>,2> _element;
  public:
    matrix2d(scalar m11=0.0,scalar m12=0.0,
	     scalar m21=0.0,scalar m22=0.0)
    {
      _element[0][0]=m11;
      _element[0][1]=m12;
      _element[1][0]=m21;
      _element[1][1]=m22;
    }
    matrix2d(const matrix2d& in)
      : _element(in._element) {}
    /*public operators*/
    //accessors
    scalar
    operator() (int i,int j) const
    {
      return _element[i][j];
    }
    scalar&
    operator() (int i, int j)
    {
      return _element[i][j];
    }
    class iterator : public niterator<std::size_t,std::size_t>
    {
      matrix2d& _m;
    public:
      iterator(noob3d::matrix2d& in,
       	       std::size_t i=0,std::size_t j=0)
	: niterator<size_t,size_t>(std::make_tuple(i,j),
				   std::make_tuple(2,2)),
	  _m(in)
      {}
      noob3d::scalar&
      operator*()
      {
	return _m(std::get<0>(t),std::get<1>(t));
      }
    };
    iterator
    begin() { return iterator(*this); };
    iterator
    end() { return iterator(*this,0,2); };

    //normal arithmetic
    matrix2d
    operator+ (const matrix2d &rval) const
    {
      matrix2d ret(*this);
      /*only unrolled for loops*/
      ret(0,0) += rval(0,0);
      ret(1,0) += rval(1,0);
      ret(0,1) += rval(0,1);
      ret(1,1) += rval(1,1);
      return ret;
    }

    matrix2d
    operator- (const matrix2d &rval) const
    {
        matrix2d ret(*this);
	/*only unrolled for loops*/
	ret(0,0) -= rval(0,0);
	ret(1,0) -= rval(1,0);
	ret(0,1) -= rval(0,1);
	ret(1,1) -= rval(1,1);
	return ret;	
    }

    matrix2d
    operator* (scalar rval) const
    {
      matrix2d ret(*this);
      /*only unrolled for loops*/
      ret(0,0) *= rval;
      ret(1,0) *= rval;
      ret(0,1) *= rval;
      ret(1,1) *= rval;
      return ret;
    }

    matrix2d
    operator/ (scalar rval) const
    {
      matrix2d ret(*this);
      /*only unrolled for loops*/
      ret(0,0) /= rval;
      ret(1,0) /= rval;
      ret(0,1) /= rval;
      ret(1,1) /= rval;
      return ret;
    }

    vector2d
    operator* (const vector2d &rval) const
    {
      vector2d ret;
      ret[0] = (*this)(0,0)*rval[0]+(*this)(0,1)*rval[1];
      ret[1] = (*this)(1,0)*rval[0]+(*this)(1,1)*rval[1];
      return ret;
    }

    matrix2d
    operator* (const matrix2d &rval) const
    {
      matrix2d ret;
      ret(0,0) = (*this)(0,0)*rval(0,0)+(*this)(0,1)*rval(1,0);
      ret(0,1) = (*this)(0,0)*rval(0,1)+(*this)(0,1)*rval(1,1);

      ret(1,0) = (*this)(1,0)*rval(0,0)+(*this)(1,1)*rval(1,0);
      ret(1,1) = (*this)(1,0)*rval(0,1)+(*this)(1,1)*rval(1,1);
      return ret;
    }

    //more efficient operators, use these

    matrix2d&
    operator= (const matrix2d &rval) //note, this is arithmetic assign!
    { _element=rval._element; return *this;}

    matrix2d&
    operator+= (const matrix2d &rval)
    {
      (*this)(0,0)+=rval(0,0);
      (*this)(1,0)+=rval(1,0);
      (*this)(0,1)+=rval(0,1);
      (*this)(1,1)+=rval(1,1);
      return *this;
    }

    matrix2d&
    operator-= (const matrix2d &rval)
    {
      (*this)(0,0)-=rval(0,0);
      (*this)(1,0)-=rval(1,0);
      (*this)(0,1)-=rval(0,1);
      (*this)(1,1)-=rval(1,1);
      return *this;
    }

    matrix2d&
    operator*= (scalar rval)
    {
      (*this)(0,0)*=rval;
      (*this)(1,0)*=rval;
      (*this)(0,1)*=rval;
      (*this)(1,1)*=rval;
      return *this;
    }
    matrix2d&
    operator*= (const matrix2d &rval)
    {
      matrix2d ret;
      ret(0,0) = (*this)(0,0)*rval(0,0)+(*this)(0,1)*rval(1,0);
      ret(0,1) = (*this)(0,0)*rval(0,1)+(*this)(0,1)*rval(1,1);

      ret(1,0) = (*this)(1,0)*rval(0,0)+(*this)(1,1)*rval(1,0);
      ret(1,1) = (*this)(1,0)*rval(0,1)+(*this)(1,1)*rval(1,1);
      *this = ret;//'tis unfortunate.
      return *this;
    }
    matrix2d&
    transpose()
    {
      std::swap((*this)(1,0),(*this)(0,1));
      return *this;
    }
  };

  class matrix3d
  {
    std::array<std::array<scalar,3>,3> _element;
  public:
    matrix3d(scalar m11=0.0,scalar m12=0.0, scalar m13=0.0,
	     scalar m21=0.0,scalar m22=0.0, scalar m23=0.0,
	     scalar m31=0.0,scalar m32=0.0, scalar m33=0.0)
    {
      _element[0][0]=m11;
      _element[0][1]=m12;
      _element[0][2]=m13;

      _element[1][0]=m21;
      _element[1][1]=m22;
      _element[1][2]=m23;

      _element[2][0]=m31;
      _element[2][1]=m32;
      _element[2][2]=m33;
    }
    matrix3d(const matrix3d& in)
      : _element(in._element) {}
    
    /*public operators*/
    //accessors
    scalar
    operator() (int i,int j) const
    {
      return _element[i][j];
    }
    scalar&
    operator() (int i, int j)
    {
      return _element[i][j];
    }
    
    class citerator : public niterator<std::size_t,std::size_t>
    {
      const matrix3d& _m;
    public:
      citerator(const noob3d::matrix3d& in,
		std::size_t i=0,std::size_t j=0)
	: niterator<size_t,size_t>(std::make_tuple(i,j),
				   std::make_tuple(3,3)),
	  _m(in)
      {}
      noob3d::scalar
      operator*()
      {
	return _m(std::get<0>(t),std::get<1>(t));
      }
      citerator
      corresponding(noob3d::matrix3d& in)
      {
	return citerator(in,std::get<0>(t),std::get<1>(t));
      }
    };
    class iterator : public niterator<std::size_t,std::size_t>
    {
      matrix3d& _m;
    public:
      iterator(noob3d::matrix3d& in,
       	       std::size_t i=0,std::size_t j=0)
	: niterator<size_t,size_t>(std::make_tuple(i,j),
				   std::make_tuple(3,3)),
	  _m(in)
      {}
      noob3d::scalar&
      operator*()
      {
	return _m(std::get<0>(t),std::get<1>(t));
      }
      iterator
      corresponding(noob3d::matrix3d& in) const
      {
	return iterator(in,std::get<0>(t),std::get<1>(t));
      }
    };    
    iterator
    begin() { return iterator(*this); };
    iterator
    end() { return iterator(*this,0,3); };
    citerator
    cbegin() { return citerator(*this); };
    citerator
    cend() { return citerator(*this,0,3); };

    //normal arithmetic
    matrix3d
    operator+ (const matrix3d &rval) const
    {
      matrix3d ret(*this);
      /*only unrolled for loops*/
      ret(0,0) += rval(0,0);
      ret(0,1) += rval(0,1);
      ret(0,2) += rval(0,2);

      ret(1,0) += rval(1,0);
      ret(1,1) += rval(1,1);
      ret(1,2) += rval(1,2);

      ret(2,0) += rval(2,0);
      ret(2,1) += rval(2,1);
      ret(2,2) += rval(2,2);
      return ret;
    }

    matrix3d
    operator- (const matrix3d &rval) const
    {
      matrix3d ret(*this);
      /*only unrolled for loops*/


      ret(0,0) -= rval(0,0);
      ret(0,1) -= rval(0,1);
      ret(0,2) -= rval(0,2);

      ret(1,0) -= rval(1,0);
      ret(1,1) -= rval(1,1);
      ret(1,2) -= rval(1,2);

      ret(2,0) -= rval(2,0);
      ret(2,1) -= rval(2,1);
      ret(2,2) -= rval(2,2);
      return ret;	
    }

    matrix3d
    operator* (scalar rval) const
    {
      matrix3d ret(*this);
      /*only unrolled for loops*/
      ret(0,0) *= rval;
      ret(0,1) *= rval;
      ret(0,2) *= rval;

      ret(1,0) *= rval;
      ret(1,1) *= rval;
      ret(1,2) *= rval;

      ret(2,0) *= rval;
      ret(2,1) *= rval;
      ret(2,2) *= rval;
      return ret;	
    }

    matrix3d
    operator/ (scalar rval) const
    {
      matrix3d ret(*this);
      /*only unrolled for loops*/
      ret(0,0) /= rval;
      ret(0,1) /= rval;
      ret(0,2) /= rval;

      ret(1,0) /= rval;
      ret(1,1) /= rval;
      ret(1,2) /= rval;

      ret(2,0) /= rval;
      ret(2,1) /= rval;
      ret(2,2) /= rval;
      return ret;	
    }

    vector3d
    operator* (const vector3d &rval) const
    {
      vector3d ret;
      ret[0] = (*this)(0,0)*rval[0]+(*this)(0,1)*rval[1]+(*this)(0,2)*rval[2];
      ret[1] = (*this)(1,0)*rval[0]+(*this)(1,1)*rval[1]+(*this)(1,2)*rval[2];
      ret[2] = (*this)(2,0)*rval[0]+(*this)(2,1)*rval[1]+(*this)(2,2)*rval[2];
      return ret;
    }

    matrix3d
    operator* (const matrix3d &rval) const
    {
      matrix3d ret;
      ret(0,0) = (*this)(0,0)*rval(0,0)+(*this)(0,1)*rval(1,0)+(*this)(0,2)*rval(2,0);
      ret(0,1) = (*this)(0,0)*rval(0,1)+(*this)(0,1)*rval(1,1)+(*this)(0,2)*rval(2,1);
      ret(0,2) = (*this)(0,0)*rval(0,2)+(*this)(0,1)*rval(1,2)+(*this)(0,2)*rval(2,2);

      ret(1,0) = (*this)(1,0)*rval(0,0)+(*this)(1,1)*rval(1,0)+(*this)(1,2)*rval(2,0);
      ret(1,1) = (*this)(1,0)*rval(0,1)+(*this)(1,1)*rval(1,1)+(*this)(1,2)*rval(2,1);
      ret(1,2) = (*this)(1,0)*rval(0,2)+(*this)(1,1)*rval(1,2)+(*this)(1,2)*rval(2,2);

      ret(2,0) = (*this)(2,0)*rval(0,0)+(*this)(2,1)*rval(1,0)+(*this)(2,2)*rval(2,0);
      ret(2,1) = (*this)(2,0)*rval(0,1)+(*this)(2,1)*rval(1,1)+(*this)(2,2)*rval(2,1);
      ret(2,2) = (*this)(2,0)*rval(0,2)+(*this)(2,1)*rval(1,2)+(*this)(2,2)*rval(2,2);
      return ret;
    }

    //more efficient operators, use these

    matrix3d&
    operator= (const matrix3d &rval) //note, this is arithmetic assign!
    { _element=rval._element; return *this;}
    matrix3d&
    operator+= (const matrix3d &rval)
    {
      (*this)(0,0)+=rval(0,0);
      (*this)(0,1)+=rval(0,1);
      (*this)(0,2)+=rval(0,2);

      (*this)(1,0)+=rval(1,0);
      (*this)(1,1)+=rval(1,1);
      (*this)(1,2)+=rval(1,2);

      (*this)(2,0)+=rval(2,0);
      (*this)(2,1)+=rval(2,1);
      (*this)(2,2)+=rval(2,2);
      return *this;
    }

    matrix3d&
    operator-= (const matrix3d &rval)
    {
      (*this)(0,0)-=rval(0,0);
      (*this)(0,1)-=rval(0,1);
      (*this)(0,2)-=rval(0,2);

      (*this)(1,0)-=rval(1,0);
      (*this)(1,1)-=rval(1,1);
      (*this)(1,2)-=rval(1,2);

      (*this)(2,0)-=rval(2,0);
      (*this)(2,1)-=rval(2,1);
      (*this)(2,2)-=rval(2,2);
      return *this;
    }

    matrix3d&
    operator*= (scalar rval)
    {
      (*this)(0,0)*=rval;
      (*this)(0,1)*=rval;
      (*this)(0,2)*=rval;

      (*this)(1,0)*=rval;
      (*this)(1,1)*=rval;
      (*this)(1,2)*=rval;

      (*this)(2,0)*=rval;
      (*this)(2,1)*=rval;
      (*this)(2,2)*=rval;
      return *this;
    }
    matrix3d&
    operator/= (scalar rval)
    {
      (*this)(0,0)/=rval;
      (*this)(0,1)/=rval;
      (*this)(0,2)/=rval;

      (*this)(1,0)/=rval;
      (*this)(1,1)/=rval;
      (*this)(1,2)/=rval;

      (*this)(2,0)/=rval;
      (*this)(2,1)/=rval;
      (*this)(2,2)/=rval;
      return *this;
    }

    matrix3d&
    operator*= (const matrix3d &rval)
    {
      matrix3d ret;
      ret(0,0) = (*this)(0,0)*rval(0,0)+(*this)(0,1)*rval(1,0)+(*this)(0,2)*rval(2,0);
      ret(0,1) = (*this)(0,0)*rval(0,1)+(*this)(0,1)*rval(1,1)+(*this)(0,2)*rval(2,1);
      ret(0,2) = (*this)(0,0)*rval(0,2)+(*this)(0,1)*rval(1,2)+(*this)(0,2)*rval(2,2);

      ret(1,0) = (*this)(1,0)*rval(0,0)+(*this)(1,1)*rval(1,0)+(*this)(1,2)*rval(2,0);
      ret(1,1) = (*this)(1,0)*rval(0,1)+(*this)(1,1)*rval(1,1)+(*this)(1,2)*rval(2,1);
      ret(1,2) = (*this)(1,0)*rval(0,2)+(*this)(1,1)*rval(1,2)+(*this)(1,2)*rval(2,2);

      ret(2,0) = (*this)(2,0)*rval(0,0)+(*this)(2,1)*rval(1,0)+(*this)(2,2)*rval(2,0);
      ret(2,1) = (*this)(2,0)*rval(0,1)+(*this)(2,1)*rval(1,1)+(*this)(2,2)*rval(2,1);
      ret(2,2) = (*this)(2,0)*rval(0,2)+(*this)(2,1)*rval(1,2)+(*this)(2,2)*rval(2,2);
      *this=ret;
      return *this;
    }
    matrix3d&
    transpose()
    {
      std::swap((*this)(1,0),(*this)(0,1));
      std::swap((*this)(2,0),(*this)(0,2));
      std::swap((*this)(1,2),(*this)(2,1));
      return *this;
    }
  };

  class matrix3d_cofactor //will generalize once I make a ndmatrix type.
  {
    const matrix3d &_m;
    std::size_t si,sj;
  public:
    matrix3d_cofactor(const matrix3d& im,std::size_t i,std::size_t j) 
      : _m(im), si(i), sj(j){}
    scalar
    operator()(std::size_t i, std::size_t j) const
    {
      if (i >= si) ++i;
      if (j >= sj) ++j;
      return _m(i,j);
    }
  };

  void
  print(const matrix3d& m)
  {
    std::cout << std::setprecision(2)
      << "|" << m(0,0) << "," << m(0,1)<< "," << m(0,2) << "|" << std::endl
      << "|" << m(1,0) << "," << m(1,1)<< "," << m(1,2) << "|" << std::endl
      << "|" << m(2,0) << "," << m(2,1)<< "," << m(2,2) << "|" << std::endl;
    
  }
  void
  print(const matrix2d& m)
  {
    std::cout << "|" << m(0,0) << "," << m(0,1)<<  "|" << std::endl
	      << "|" << m(1,0) << "," << m(1,1)<<  "|" << std::endl;
  }
  void
  print(const matrix3d_cofactor& m)
  {
    std::cout << "|" << m(0,0) << "," << m(0,1)<<  "|" << std::endl
	      << "|" << m(1,0) << "," << m(1,1)<<  "|" << std::endl;
  }
  void
  print(matrix3d::iterator i)
  {
    std::cout << "(" << std::get<0>(i.t) << "," << std::get<1>(i.t) << ")" << std::endl;
  }

  inline matrix3d_cofactor
  cofactor(const matrix3d& m, matrix3d::iterator i)
  {
    return matrix3d_cofactor(m,std::get<0>(i.t),std::get<1>(i.t));
  }

  inline matrix3d_cofactor
  cofactor(const matrix3d& m, matrix3d::citerator i)
  {
    return matrix3d_cofactor(m,std::get<0>(i.t),std::get<1>(i.t));
  }
  template <typename Matrix2d>
  scalar
  det2d(const Matrix2d& m)
  {
    return m(0,0)*m(1,1)-m(1,0)*m(0,1);
  }
  
  inline scalar
  det(const matrix2d& m)
  {
    return det2d<matrix2d>(m);
  }
  inline scalar
  det(const matrix3d_cofactor& m)
  {
    return det2d<matrix3d_cofactor>(m);
  }
  
  inline scalar
  det(const matrix3d& m)
  {
    return 
      m(0,0)*m(1,1)*m(2,2)
      +m(0,1)*m(1,2)*m(2,0)
      +m(0,2)*m(1,0)*m(2,1)
      -m(0,0)*m(1,2)*m(2,1)
      -m(0,1)*m(1,0)*m(2,2)
      -m(0,2)*m(1,1)*m(2,0);
  }
  inline matrix3d
  inverse(const matrix3d& m)
  {
    matrix3d o;
    int sign=1;
    for(auto i = o.begin(); i != o.end(); ++i)
      {
	*i = sign*det(cofactor(m,i));
	sign*=-1;
      }
    return o.transpose()/det(m);
  }
  inline matrix2d
  inverse(const matrix2d& m)
  {
    return matrix2d( m(1,1),-m(0,1),
		    -m(1,0), m(0,0))/det(m);
  }
} //namespace noob3d
#endif //_NOOB3D_MATRIX_H_

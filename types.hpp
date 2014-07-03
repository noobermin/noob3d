#ifndef _NOOB3D_TYPES_H_
#define _NOOB3D_TYPES_H_

#include <vector>
#include <tuple>
#include <exception>
#include <stdexcept>
#include <type_traits>
#include <iterator>
namespace noob3d
{
  //an iterator for iterating over a container
  //overloading the (i,j,...) operator
  //with types Indices..., naturally.
  //starts from i and increases.
  template <typename...Indices>
  struct niterator
  {
    //increment
    //tail end
    template <std::size_t I=0> 
    typename std::enable_if<I == sizeof...(Indices)-1, niterator<Indices...>&>::type 
    operator++()
    {
      ++std::get<I>(t);
      return *this;
    }
    template <std::size_t I=0>
    typename std::enable_if< I < sizeof...(Indices)-1, niterator<Indices...>&>::type 
    operator++()
    {
      ++std::get<I>(t);
      if (std::get<I>(t) == std::get<I>(s))
	{
	  std::get<I>(t)=0;
	  operator++<I+1>();
	}
      return *this;
    }
    //comparator
    //tail end
    template <std::size_t I=0> 
    typename std::enable_if<I == sizeof...(Indices), bool>::type 
    operator==(niterator<Indices...> in) const
    {
      return true;
    }
    template <std::size_t I=0>
    typename std::enable_if< I < sizeof...(Indices), bool>::type 
    operator==(niterator<Indices...> in) const
    {
      if (std::get<I>(t) == std::get<I>(in.t))
	return operator==<I+1>(in);
      else
	return false;
    }

    std::tuple<Indices...> t,s;
    niterator(std::tuple<Indices...> indices,
	      std::tuple<Indices...> size)
      : t(indices), s(size) {}
    niterator(const niterator<Indices...>& in)
      : t(in.t), s(in.s) {}
    niterator
    operator++(int)
    {niterator<Indices...> ret(*this); this->operator++(); return ret;}
    bool operator!=(const niterator<Indices...>& in) const
    { return !operator==(in);}
  };
  //now, for my secrect sauce :)
  template <typename MultiIterator, typename... Indices>
  MultiIterator
  up_from(const MultiIterator& i, Indices... indices)
  {
    MultiIterator ret(i);
    ret.up(indices... );
    return ret;
  }
  
  template <typename MultiIterator, typename... Indices>
  MultiIterator
  down_from(const MultiIterator& i, Indices... indices)
  {
    MultiIterator ret(i);
    ret.up(indices... );
    return ret;
  }
}//namespace noob3d

#endif //_NOOB3D_TYPES_H_

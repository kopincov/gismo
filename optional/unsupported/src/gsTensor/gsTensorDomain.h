//
// Purpose: iterates over elements, samples domain, 
// produces a quadrature rule on the domain
//
#pragma once

#include <ostream>

#include <gsCore/gsDomain.h>

namespace gismo
{

  /** 
      @brief Class for iterating elements,
      sample domain and produces a quadrature rule on the domain
  */

  
template<class T, unsigned d>
class gsTensorDomain : public gsDomain<T>
{

public:

  /// Default empty constructor
  gsTensorDomain() { };

  ~gsTensorDomain() { }; //destructor


public:



// Data members
private:


}; // TensorDomain gsClass


//////////////////////////////////////////////////
//////////////////////////////////////////////////




}; // namespace gismo

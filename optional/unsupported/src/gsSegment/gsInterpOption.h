/** @file gsInterpOption.h

    @brief Provides declaration of gsInterpOption class. Initialization for
    default options and enabling manual adjustments for tensor product spline
    surface interpolation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D.-M. Nguyen
*/

#pragma once

#include <iostream>

namespace gismo
{
template<class T>
class gsInterpOption
{
  
// Data members
public:
  T wReg; // weight: regularization
  T wdPoint; // weight: discrete data: interpolating points
  T wdNormal; // weight: discrete data: normals
  int ndInPoint; // number: discrete
  T wcPoint; // weight: continuous data: interpolating points
  T wcNormal; // weight: continuous data: normals
  
public:

  /// Default empty constructor
  gsInterpOption(){ }

  /// Constructing using a spline degree
  gsInterpOption(const int & splineDeg)
  {
      wReg = 1;
      wdPoint = 1;
      wdNormal = 1;
      ndInPoint = splineDeg-1;      
      wcPoint = 1;
      wcNormal = 1;
  }

  /// Destructor
  virtual ~gsInterpOption() {}
  
   void update(std::string const & field,T const & value)
   {
     if (field=="wReg") wReg = value;
     if (field=="wdPoint") wdPoint = value;
     if (field=="wdNormal") wdNormal = value;
//      else
//      {
//        std::cout << std::endl << "The intered name does not match any of the members of the class gsInterpOption" << std::endl;
//        exit(1);
//      };

   }
  
  /// Prints the object as a string
  ///
  /// \param os The output stream.
  virtual std::ostream &print(std::ostream &os) const 
  {
    os<<"gsInterpOption:\n";
    //os<<"InfRep: "<< InfRep << std::endl;
    return os;       
  }
  


};

} // namespace


//////////////////////////////////////////////////
//////////////////////////////////////////////////
/// Define the functionality of the operator << acting on a class
// template<class T>
// std::ostream &operator<<(std::ostream &os, const gsInterpOption<T>& b)
// {return b.print(os); };

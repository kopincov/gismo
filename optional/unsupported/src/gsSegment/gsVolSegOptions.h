/** @file gsVolSegOptions.h

    @brief Provides declaration of gsVolSegOptions class. Options for volume
    segmentation of a solid given by its trimmed surfaces (gsSolid).
    loop.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D.-M. Nguyen
*/

#pragma once

#include <gsUtils/gsGraph/gsGraphOptions.h>

namespace gismo
{
  
template<class T>
class gsVolSegOptions: gsGraphOptions<T>
{
  
// Data members
public:
  T costConvexEdge; // cost of convex edges
  T costNonConvexEdge; // cost of non-convex edges
  T costAuxEdge; // penalty cost of auxiliary edges
  bool costAuxEdgeAngle; // if true then costs of auxiliary edges = penalty cost + angle discrepancy
  bool DiffFaces; // if true then edges in a path are on different faces
  bool PlanarityCost; // if true then add planarity value to the cutting algorithm  
  int PlanarityCostType; // = 1: sum of unsigned dis/objectDimemsion*100
  T PlanarityCostWeight;
  
public:

  /// Default empty constructor
  gsVolSegOptions(): gsGraphOptions<T>() 
  { 
      costConvexEdge = 1;
      costNonConvexEdge = 0;
      costAuxEdge = 10;
      costAuxEdgeAngle = true;
      DiffFaces = true;
      PlanarityCost = true;
      PlanarityCostType = 1;
      PlanarityCostWeight = 1;
  }

  /// Destructor
  ~gsVolSegOptions() {}
  
  /// Set an option.
  /// \param field Name of the option to set. Can be: "costConvexEdge",
  /// "costNonConvexEdge", "costAuxEdge", "costAuxEdgeAngle", "DiffFaces",
  /// "PlanarityCost", "PlanarityCostWeight", "PlanarityCostType".
  /// \param value Value for this option.
  void update(std::string field,T value)
  {
    if (field=="costConvexEdge")
      costConvexEdge = value;
    else if (field=="costNonConvexEdge")
      costNonConvexEdge = value;    
    else if (field=="costAuxEdge")
      costAuxEdge = value;
    else if (field=="costAuxEdgeAngle")
      costAuxEdgeAngle = value;
    else if (field=="DiffFaces")
      DiffFaces = value;
    else if (field=="PlanarityCost")
      PlanarityCost = value;
    else if (field=="PlanarityCostWeight")
      PlanarityCostWeight = value;
    else if (field=="PlanarityCostType")
      PlanarityCostType = 1;        
    else
    {
      std::cout << std::endl << "The intered name does not match any of the members of the class gsVolSegOptions" << std::endl;
      exit(1);
    };
      
  }
  

  /// Prints the object as a string
  /// \param os Output stream to print to.
  virtual std::ostream &print(std::ostream &os) const 
  {
    os<<"gsVolSegOptions:\n";
    os<<"costConvexEdge: "<< costConvexEdge << std::endl;
    os<<"costAuxEdge: "<< costAuxEdge << std::endl;
    os<<"DiffFaces: "<< DiffFaces << std::endl;
    return os;       
  };  
  


}; 


//////////////////////////////////////////////////
//////////////////////////////////////////////////
/// Define the functionality of the operator << acting on a class
template<class T>
std::ostream &operator<<(std::ostream &os, const gsVolSegOptions<T>& b)
{return b.print(os); };




}; // namespace gismo

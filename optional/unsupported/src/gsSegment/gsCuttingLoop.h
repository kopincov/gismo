/** @file gsCuttingLoop.h

    @brief Provides declaration of gsCuttingLoop class. A loop of curves used
    for determining a cutting surface which segments the volume into two.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D.-M. Nguyen
*/

#pragma once

#include <iostream>
#include <gsUtils/gsGraph/gsGraphPath.h>
#include <gsModeling/gsSolid.h>

namespace gismo
{

template <class T > class gsSolid;
template<class T>
class gsCuttingLoop: public gsGraphPath<T>
{
private:
  gsSolid<T>* solid;
  
public:

  /// Default empty constructor
  gsCuttingLoop(): gsGraphPath<T>() { };
  
  gsCuttingLoop(gsGraphPath<T>* path, gsSolid<T>* solidI): gsGraphPath<T>(path->getSource(),path->getTarget(),path->getPrevMat(),path->getLength())
  {
    solid = solidI;
  };  
  
  
  ~gsCuttingLoop() { }; //destructor


public:
  /// Return number of vertices
  int size() const
  {
    std::vector<unsigned> path=this->computePath();
    return path.size();
  };
  
  /// Return coordinates of the vertices of the cutting loop
  gsMatrix<T> getVertices() const;
  
  /// sample points on the cutting loop, output: a pointer to a matrix of size 3 x points*numVertices
  /// There MUST be half edges created already for auxiliary edges
  ///
  /// \param npoints Number of points to sample.
  /// \param numEndPoints Include this many endpoints in the sample (0, 1 or 2).
  gsMatrix<T> sample(int npoints, int numEndPoints=2) const;
  
  /// sample normals on the cutting loop wrt the order of vertices in the cutting loop
  /// this order depends the assigment of target vertex and source vertex when finding a cutting loop
  /// There MUST be half edges created already for auxiliary edges
  ///
  /// \param npoints Number of points to sample.
  gsMatrix<T> sampleNormal(int npoints) const;


}; // class gsCuttingLoop


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCuttingLoop.hpp)
#endif



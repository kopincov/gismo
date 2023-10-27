/** @file gsCuttingLoop.hpp

    @brief Provides implementation of gsCuttingLoop class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D.-M. Nguyen
*/

#pragma once

#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsModelingUtils.hpp>

namespace gismo
{

template <class T>
gsMatrix<T> gsCuttingLoop<T>::getVertices() const
{
  std::vector<unsigned> path=this->computePath();
  size_t nVert = path.size();

  std::rotate( path.begin(), path.end()-1, path.end() );

  gsMatrix<T> vert(3,nVert);
  for (size_t i=0;i!=nVert;i++)
  {
    vert.col(i) = solid->getVertexFromID(path[i])->getCoordinate();
    //vert.col(nVert-1-i) = solid->getVertexFromID(path[i])->getCoordinate();
  };
  return vert;
};

template <class T>
gsMatrix<T> gsCuttingLoop<T>::sample(int npoints, int numEndPoints) const
{
  assert(npoints>=2);
  assert(numEndPoints>=0 && numEndPoints<=2);     
  bool const sameSpeed=true; // assume that the two trimming curves have the same speed (but in different direction)
  int consideredCurve=1; // =1: average of two trimming curves
  assert(sameSpeed==true && consideredCurve==1); // if not, TODO: reparametrize the two curves to have constant speed curves
  int np; // new number of points
  switch (numEndPoints)
  {
    case (0):
      np = npoints-2;
      break;
    case (1):
      np = npoints-1;
      break;	
    case (2):
      np = npoints;
  default:
      np = npoints;
  }
  
  std::vector<unsigned> path=this->computePath();
  int nVert = path.size();  

  std::rotate( path.begin(), path.end()-1, path.end() );
      
  gsMatrix<T> image(3,nVert*np);
  unsigned v1,v2;
  
  int ind=0;
  gsSolidHalfEdge<T> * he;
  int tloopInd,tloopIndMate;
  gsMatrix<T> uCols, uColsMate;
  
  if (sameSpeed==true && consideredCurve==1)
  {
    for (int i=0;i<=nVert-1;i++)
    {
      v1 = path[i];
      v2 = path[(i+1) % nVert];
      he = solid->getVertexFromID(v1)->getHalfEdge(solid->getVertexFromID(v2)); // he must exist 
      tloopInd = he->trimLoopInd(he->eps * 1e+8);
      tloopIndMate = he->mate->trimLoopInd(he->eps * 1e+8);
      he->face->surf->sampleCurve_into(he->loopN(), tloopInd, npoints, uCols);
      he->mate->face->surf->sampleCurve_into(he->mate->loopN(), tloopIndMate, npoints, uColsMate);
	  
      image.middleCols( ind * np,np ) = .5*uCols + .5*flipLR<T>(uColsMate);
      ind++;
    }
  }

  return image;
}


template <class T>
gsMatrix<T> gsCuttingLoop<T>::sampleNormal(int npoints) const
{
  assert(npoints>=2);
  bool const sameSpeed=true; // assume that the two trimming curves have the same speed (but in different direction)
  int consideredCurve=1; // =1: average of two trimming curves
  assert(sameSpeed==true && consideredCurve==1); // if not, TODO: reparametrize the two curves to have constant speed curves
  int np = npoints;      
  
  std::vector<unsigned> path=this->computePath();
  int nVert = path.size();

  std::rotate( path.begin(), path.end()-1, path.end() );

  gsMatrix<T> image(3,nVert*np);
  image.setZero();
  unsigned v1,v2;
  
  int ind=0;
  gsSolidHalfEdge<T> * he;
  int tloopInd,tloopIndMate;
  gsMatrix<T> uCols, uColsMate, pointingIn; // pointingIn: vectors pointing into the "center" of the cutting surface
  
  if (sameSpeed==true && consideredCurve==1)
  {
    for (int i=0;i<=nVert-1;i++)
    {
      v1 = path[i];
      v2 = path[(i+1) % nVert];
      he = solid->getVertexFromID(v1)->getHalfEdge(solid->getVertexFromID(v2)); // he must exist 
      tloopInd = he->trimLoopInd(he->eps * 100);
      tloopIndMate = he->mate->trimLoopInd(he->eps * 100);
      uCols = he->face->surf->sampleNormal(he->loopN(),tloopInd,npoints);
      uColsMate = he->mate->face->surf->sampleNormal(he->mate->loopN(),tloopIndMate,npoints);
      uColsMate = flipLR<T>(uColsMate);
      pointingIn = - (uCols + uColsMate);
      // get tangent vectors of the trimming curve at a list of uniformly spaced points 
      gsMatrix<T> tangents = he->face->surf->trimCurTangents(he->loopN(),tloopInd,npoints);
      //gsMatrix<T> tet = crossNorm2Mat<T>(tangents,pointingIn);
      image.middleCols( ind * np,np ) = crossNorm2Mat<T>(tangents,pointingIn);
      ind++;
    };
  };  
  
  return image;
};  


}; // namespace gismo

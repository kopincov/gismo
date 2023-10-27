/** @file gsVolumeSegment.h

    @brief Provides declaration of gsVolumeSegment class. Provides
    functionality for the segmentation of a boundary represented solid.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D.-M. Nguyen, M. Pauley
*/

//gsinterpolate// TODO

// CUTTING SURFACE
// - take care of the fact that two trimming curves representing one solid edge should have the same speed
// - take care of scaling control points of trimming curves in parameter domains to between 0 and 1
// - improve the code: operator= for Bsplines, multiplication of matrices, safe evaluation, columewise multiplication

// CHANGES DUE TO THE VOLUME SPLIT
// - change edges of a segmenting loop from nonconvex to convex once the segmentation is done.
// - update adjacent matrix, especially for new vertices, and new faces
// - update input for shortest path, as two vertices does not define a halfedge,
// as there are several HEs in different sub-volumes connecting the two vertices

// COST OF AUX EDGES
// - improve cost for aux edge: 
// + check if it emanates from a non-convex vertex,
// + check if it split the face into 4-sided face
// + define a tolerance on the magnitude of a tentatively cut angle

// PRACTICAL
// - check with Angelos if it is possible to have suggesting dropping list when type *make* in build folder
// - check with Angelos how stable the computations of derivatives at end points are

// IMPLIMENTATION: 

// NOTE:

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsModeling/gsSolid.h>
#include <gsSegment/gsVolSegOptions.h>
#include <gsUtils/gsGraph/gsGraph.h>
#include <gsSegment/gsCuttingLoop.h>
#include <gsSegment/gsInterpOption.h>

namespace gismo
{

  /** 
      @brief An algorithm for volume segmentation of a solid given by its trimmed surfaces (gsSolid). 
  */

  
template<class T>
class gsVolumeSegment: public gsGraph<T>
{
public:
  typedef typename gsGraph<T>::gsMatrixT gsMatrixT;

public:

  /// Default empty constructor
  gsVolumeSegment(): gsGraph<T>() {}
  
  /// Construct the object for a given gsSolid *solid* and options *opt*
  gsVolumeSegment(gsSolid<T> * sl); 
  
  /// For a given halfedge \a he, find a cutting loop from the source of *he* to its target
  /// \param he The half-edge whose ends are used as endpoints for the search.
  gsCuttingLoop<T>* FindCuttingLoop(gsSolidHalfEdge<T> * he)
  {
    unsigned vert1 = he->source->getId();
    unsigned vert2 = he->target()->getId();
    gsGraphPath<T>* path = this->Dijkstra(vert1,vert2,he);
    gsCuttingLoop<T> * cloop = new gsCuttingLoop<T>(path,this->solid);
    delete path;
    return cloop;
  }; 
  
  /// For given vertices with IDS: \a vertexID1 and \a vertexID2, find a cutting loop from vertexID1 to vertexID2
  /// \param vertexID1 Start point for the search.
  /// \param vertexID2 End point for the search.
  gsCuttingLoop<T>* FindCuttingLoop(int vertexID1, int vertexID2)
  {
    gsSolidHalfEdge<T>* he = this->solid->getVertexFromID(vertexID1)->getHalfEdge( this->solid->getVertexFromID(vertexID2) );
    return this->FindCuttingLoop( he );
  }; 
  
public:


  /// Intepolating a cutting loop into a trimmed surface with the curves being its trimming curves
  /// \param cloop The cutting loop that will be interpolated to construct the surface.
  /// \param curveLoop Trimming loop for the domain.
  /// \param kv1 Knot vector in the first dimension.
  /// \param kv2 Knot vector in the second dimension.
  /// \param intpOpt a gsInterpOption object providing the options for the interpolation
  /// \param force_normal Used for a special case. Normally should be false.
  gsTrimSurface<T>* cuttingSurface(  gsCuttingLoop<T> const & cloop, const gsCurveLoop<T> & curveLoop, gsKnotVector<T> const & kv1,
                    gsKnotVector<T> const & kv2,gsInterpOption<T> & intpOpt, bool force_normal=false) const;

  /// Intepolating descrete and smooth data by a tensor-product spline surface
  /// TODO: will be disscussed if it should be moved to some "interpolating place"
  typedef gsMatrix<T> pointEval(gsMatrix<T> preImage);
  typedef gsMatrix<T> normalEval(gsMatrix<T> preImage);

  /// Prints the object as a string.
  /// \param ost The output stream to print to.
  std::ostream &print(std::ostream &ost) const
      {
	  ost << "gsVolumeSegment: \n";	  	  
	  ost << "Adjacent matrix: \n" <<  this->adjMatrix << std::endl;
	  ost << this->vsOption;
	  
	  return ost;
      };    
      
  /// Cost of auxiliary edge
  T costAuxiliaryEdge(T angle,T angle1,T & angle2, T angleT,T angle1T,T & angle2T) const;
  
  /// Create an appropriate trimming loop from a cutting loop. Assumes that the
  /// edges already exist - in particular, faces have been split to make any
  /// auxiliary edges into real edges.
  /// \param solid The solid that is being segmented.
  /// \param source The cutting loop to use to construct the trimming loop.
  static gsCurveLoop<T> *curveLoopFromCuttingLoop(const gsSolid<T> & solid, const gsCuttingLoop<T> & source);
  
  /// choose an appropriate convex polygon that this surface can be mapped to.
  /// pass regular=false to choose angles based on the surface, or regular=true
  /// to just generate the vertices of a regular polygon.
  /// \param trimSurface The trimmed surface whose domain is used as the domain
  /// of the map.
  /// \param regular Set to true to use a regular polygon as the image.
  /// \param[out] Computed corners of the polygon.
  static void chooseConvexPolygon(gsTrimSurface<T> &trimSurface, bool regular, gsMatrix<T> &corners);

  /// trace a curve, using either mean value interpolation or harmonic functions. This function
  /// tries several points along a straight line to use as starting points. It keeps trying
  /// until none of the results are NaN, or it runs out of points to try. (In the latter case it
  /// fails an assertion.)
  /// \param trimSurface The trimmed surface whose domain is used as the domain
  /// of the map.
  /// \param cornerPoints The corners of the polygonal image.
  /// \param templ A template containing the curve to trace.
  /// \param firstTrialPoint A 2d vector - an initial point for the search.
  /// \param step A 2d vector - displacement of each sample point from the previous.
  /// \param maxTrials Maximum number of trials before giving in on convergence.
  /// \param harmonic Set true to use harmonic maps, false to use mean value.
  /// \param n_fit_points Number of sample points.
  /// \param traceTolerance Tolerance to use for convergence.
  static void attemptTraceCurve(gsTrimSurface<T> &trimSurface, const gsMatrix<T> cornerPoints,
                                gsTemplate<T> &templ, const gsMatrix<T> &firstTrialPoint, const gsMatrix<T> &step,
                                int maxTrials, bool harmonic, int n_fit_points, T traceTolerance, gsMatrix<T> &result);

  /// find a curve to split a single trimSurface in two
  /// \param trimSurface The trimmed surface to split.
  /// \param idxStart The index of the corner to use as the start point.
  /// \param idxEnd The index of the corner to use as the end point.
  /// \param w_reg Weighting for regularity of the curve.
  /// \param w_app Weighting for accuracy of the curve's approximation of the
  /// sample points
  /// \param n_fit_points Number of points to sample.
  /// \param traceTolerance Tolerance to use for convergence when tracing the
  /// curve.
  /// \param harmonic Set true to use harmonic maps, false to use mean value.
  /// \param regular Set to true to use a regular polygon as the image.
  /// \param curveFraction Middle fraction of the curve to use for fitting.
  /// \param interpDegree Degree of the new spline curve.
  /// \param imgCurve null or a curve in the image to trace.
  /// \param[out] cuttingCurve The resulting new curve.
  /// \param[out] fitPoints The traced points that the curve is fit to.
  static void chooseCuttingCurve(gsTrimSurface<T> &trimSurface, index_t idxStart, index_t idxEnd, T w_reg, T w_app,
                                 int n_fit_points, T traceTolerance, bool harmonic, bool regular, T curveFraction, int interpDegree,
                                 gsBSpline<T> **imgCurve, gsBSpline<T> &cuttingCurve, gsMatrix<T> &fitPoints);

  /// Split faces as necessary so that all edges in \a cuttingLoop will
  /// actually be edges in the solid. Modifies \a sl in-place. When this
  /// function needs to split a face, it maps it to a convex polygon via
  /// either transfinite mean value interpolation or harmonic functions.
  /// \param sl solid to split
  /// \param cuttingLoop a gsCuttingLoop whose edges we will trace round, splitting all necessary faces
  /// \param w_reg weighting for curve regularity
  /// \param w_app weighting for approximation of the traced curve
  /// \param n_fit_points number of points to fit when tracing the curve
  /// \param traceTolerance tolerance for using Newton's method when tracing the curve
  /// \param harmonic If true, use harmonic functions. If false, use mean value interpolation.
  /// \param curveFraction what fraction of the curve to trace (use a number < 1 to avoid getting too close to the boundary)
  static void splitFaces(gsSolid<T> & sl, const gsCuttingLoop<T> &cuttingLoop, T w_reg, T w_app,
                         int n_fit_points, T traceTolerance, bool harmonic, T curveFraction);

}; // class gsVolumeSegment

}; // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsVolumeSegment.hpp)
#endif

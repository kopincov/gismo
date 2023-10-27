/** @file gsVolumeSegment.hpp

    @brief Provides implementation of gsVolumeSegment class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): D.-M. Nguyen, M. Pauley
*/

#pragma once

#include <gsModeling/gsTemplate.h>
#include <gsSegment/gsMVInterpolation.h>
#include <gsCore/gsMath.h>

#include <gsModeling/gsTraceCurve.hpp>
#include <gsModeling/gsModelingUtils.hpp>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsBSpline.h>

namespace gismo {  
  
template<class T>
gsVolumeSegment<T>::gsVolumeSegment(gsSolid<T> * sl): gsGraph<T>()
{ 
  this->solid = sl;
  gsMatrixT  adjM(sl->numVertices,sl->numVertices); // *adjM is a zero matrix by default
  adjM.setOnes();
  adjM = - adjM;
  std::vector< gsSolidHeVertex<T>* > vt;
  T angle, angle1, angle2, angleT, angle1T, angle2T;
  T cost;  
  for (typename std::vector< gsSolidHalfFace<T> * >::const_iterator fp=sl->face.begin(); fp!=sl->face.end(); ++fp)
  {
    // Cost of existing edges
    std::vector< gsSolidHalfEdge<T>* > he = (*fp)->getHalfEdges();
    for (unsigned i=0; i<=he.size()-1;i++)
    {
      if ( he[i]->is_convex )
        adjM(he[i]->source->getId(),he[i]->next->source->getId()) = this->volSegOp().costConvexEdge;
      else
        adjM(he[i]->source->getId(),he[i]->next->source->getId()) = this->volSegOp().costNonConvexEdge;
    }
    
    // Cost of auxiliary edges
    vt.clear();
    vt = (*fp)->getVertices();
    // For each pair of vertices
    for (unsigned i=0; i<=(vt.size()-1)-2; i++)
    {
        for (unsigned j=i+2; j<=(vt.size()-1); j++)
        {
            if (vt[i]->getId()==2 && vt[j]->getId()==13)
                std::cout<< " stop here";
            // if i and j are NOT adjacent to each other
            if ( (i==0 && j==vt.size()-1)==false )
            {
                adjM(vt[i]->getId(),vt[j]->getId()) = this->volSegOp().costAuxEdge;
                adjM(vt[j]->getId(),vt[i]->getId()) = this->volSegOp().costAuxEdge;
                if (this->volSegOp().costAuxEdgeAngle)
                {
                    (*fp)->surf->cuttingAngles(i,j,&angle,&angle1,&angle2,this->solid->isFaceCCW());
                    (*fp)->surf->cuttingAngles(j,i,&angleT,&angle1T,&angle2T,this->solid->isFaceCCW());
                    cost = costAuxiliaryEdge(angle,angle1,angle2,angleT,angle1T,angle2T);
                    adjM(vt[i]->getId(),vt[j]->getId()) += cost;
                    adjM(vt[j]->getId(),vt[i]->getId()) += cost;
                };
            };
        };

    };
  }
  this->adjMatrix.swap( adjM );
}

/// Define the functionality of the operator << acting on a class
template<class T>
std::ostream &operator<<(std::ostream &os, const gsVolumeSegment<T>& b)
{return b.print(os); };

template<class T>
T gsVolumeSegment<T>::costAuxiliaryEdge(T angle,T angle1,T & angle2, T angleT,T angle1T,T & angle2T) const
{
  T cost,costT;
  //let call V1, and V2 be the vectors defining the *angle*, V is the vector tangent to the line connecting two opposite vertices on a face (in parameter domain)
  //Viewing from the direction where the orientation of the face is CCW, if V is on the right side of the line bisecting the *angle* then we say V is close to V2
  //Viewing from the direction where the orientation of the face is CCW, check if V is in between V1 and V2
  bool isVcloseToV2=true;
  bool isVcloseToV2T=false;
  if ( (angle1<=angle) && (angle2<=angle) ) cost = math::abs(angle1-angle2)/angle;
  else 
  {
    isVcloseToV2 = angle2 + angle/2 <= EIGEN_PI;
    if (isVcloseToV2==false) angle2 = 2*EIGEN_PI-angle2;
    cost = (angle + 2*angle2)/angle;
  };
  if ( (angle1T<=angleT) && (angle2T<=angleT) ) costT = math::abs(angle1T-angle2T)/angleT;
  else
  {
    isVcloseToV2T = angle2T + angleT/2 <= EIGEN_PI;    
    if (isVcloseToV2T==false) angle2T = 2*EIGEN_PI-angle2T;
    costT = (angleT + 2*angle2T)/angleT;
  };
  // if both V and VT is close to V2 (then the tentatively cutting line is on the two opposite sides of the two cut angles), then we discorage this line
  // the cost will account for the angle between V2 and the line bisecting the angle
  if ( ( (angle1<=angle) && (angle2<=angle) == false ) && ( (angle1T<=angleT) && (angle2T<=angleT) == false ) & (isVcloseToV2 != isVcloseToV2T) )
  {
    cost += EIGEN_PI - (angle2 + angle/2)/angle;
    costT += EIGEN_PI - (angle2T + angleT/2)/angleT;
  };
  return cost + costT;                
}

template<class T>
gsCurveLoop<T> * gsVolumeSegment<T>::curveLoopFromCuttingLoop(const gsSolid<T> & solid, const gsCuttingLoop<T> & source)
{
  // Look up the vertices, edges and faces.
  std::vector<unsigned> vertIdxs = source.computePath(source.getPrevMat(), source.getTarget());
  // rearrange curve so it is consistent with the output of cuttingLoop::getVertices()
  std::rotate( vertIdxs.begin(), vertIdxs.end()-1, vertIdxs.end() );

  int n = vertIdxs.size(); // number of vertices in the cutting loop
  std::vector<gsSolidHeVertex<T> *> verts;
  std::vector<gsSolidHalfEdge<T> *> edgesA;
  std::vector<gsSolidHalfEdge<T> *> edgesB;
  for(int i = 0; i < n; i++)
  {
    verts.push_back(solid.vertex[vertIdxs[i]]);
  }
  for(int i = 0; i < n; i++)
  {
    edgesA.push_back(verts[i]->getHalfEdge(verts[(i + 1) % n]));
    edgesB.push_back(edgesA[i]->mate);
  }
  // Compute the turning angles at each vertex, using the tangent vectors
  // of the splines. We assume that the intersection between two faces
  // is approximately the average of the two trimming curves.
  std::vector<T> angles;
  std::vector<bool> isConvex;
  for(int i = 0; i < n; i++)
  {
    // get average (non-unit) tangent vector of the two splines,
    // and use it to compute the angle. We compute
    //   edgeOutward0 =  (normalized sum of the two faces meeting on the
    //                    previous edge)
    //   edgeOutward1 =  (same for the next edge)
    //   vertexOutward = edgeOutward0 + edgeOutward1.
    // We use this to figure out the correct sign for the angle.
    // TODO: each step of the loop, we repeat some calculations from
    // the previous step, so there are some savings we can make.
    gsSolidHalfFace<T> *f0A = edgesA[(i + n - 1) % n]->face;
    gsSolidHalfFace<T> *f0B = edgesB[(i + n - 1) % n]->face;
    gsSolidHalfFace<T> *f1A = edgesA[i]->face;
    gsSolidHalfFace<T> *f1B = edgesB[i]->face;
    
    int vertIdx0A = f0A->indexOfVertex(verts[i]);
    int vertIdx0B = f0B->indexOfVertex(verts[i]);
    int vertIdx1A = f1A->indexOfVertex(verts[i]);
    int vertIdx1B = f1B->indexOfVertex(verts[i]);

    gsVector3d<T> vel0A = - f0A->surf->derivatives(vertIdx0A) * f0A->surf->TangentCoefs_prev(vertIdx0A);
    gsVector3d<T> vel0B = - f0B->surf->derivatives(vertIdx0B) * f0B->surf->TangentCoefs_next(vertIdx0B);
    gsVector3d<T> vel1A =   f1A->surf->derivatives(vertIdx1A) * f1A->surf->TangentCoefs_next(vertIdx1A);
    gsVector3d<T> vel1B =   f1B->surf->derivatives(vertIdx1B) * f1B->surf->TangentCoefs_prev(vertIdx1B);

    // Cosine formula.
    gsVector3d<T> vel0 = vel0A + vel0B;
    gsVector3d<T> vel1 = vel1A + vel1B;
    T cosAngle = vel1.dot(vel0) / vel0.norm() / vel1.norm();
    // adjust for numerical error that could cause angle to be nan
    if(cosAngle <= T(-1))
    {
      vel0A = - f0A->surf->derivatives(vertIdx0A) * f0A->surf->TangentCoefs_prev(vertIdx0A);
      vel0B = - f0B->surf->derivatives(vertIdx0B) * f0B->surf->TangentCoefs_next(vertIdx0B);
      vel1A =   f1A->surf->derivatives(vertIdx1A) * f1A->surf->TangentCoefs_next(vertIdx1A);
      vel1B =   f1B->surf->derivatives(vertIdx1B) * f1B->surf->TangentCoefs_prev(vertIdx1B);
    }
    if(cosAngle > T(1)) cosAngle = T(1);
    using std::acos;
    T angle = acos(cosAngle);
    GISMO_ASSERT(! math::isnan(angle), "Angle is not a number");
    angles.push_back(angle);
    
    gsVector3d<T> edgeOutward0  = (f0A->surf->cornerNormal(vertIdx0A) + f0B->surf->cornerNormal(vertIdx0B)).normalized();
    gsVector3d<T> edgeOutward1 = (f1A->surf->cornerNormal(vertIdx1A) + f1B->surf->cornerNormal(vertIdx1B)).normalized();
    gsVector3d<T> vertexOutward = edgeOutward0 + edgeOutward1;

    T dotProd = (vel1.normalized() - vel0.normalized()).dot(vertexOutward);
    if(dotProd < 0)
    {
      isConvex.push_back(true);
    }
    else
    {
      isConvex.push_back(false);
    }

  }
  
  // Compute the distances in 3d
  // TODO: find lengths of splines in 3D (currently we are just using the
  // distance between the vertices)
  std::vector<T> lengths;
  for(index_t i = 0; i < n; i++)
  {
    lengths.push_back((verts[(i + 1) % n]->coords - verts[i]->coords).norm());
  }

  return new gsCurveLoop<T>(angles, lengths, isConvex, T(0.1), true);
}


template<class T>
void gsVolumeSegment<T>::chooseConvexPolygon(gsTrimSurface<T> &trimSurface, bool regular, gsMatrix<T> &corners)
{
  // calculate angles, cap below at zero. calculate a good choice for lengths
  size_t n = trimSurface.domain().outer().size();
  std::vector<T> angles, lengths;
  for(size_t i = 0; i < n; i++)
  {
    if(regular)
    {
      angles.push_back(T(2) * EIGEN_PI / n);
      lengths.push_back(1);
    }
    else
    {
      // get final velocity of the curve before this angle
      gsMatrix<T> velPrev = -trimSurface.TangentCoefs_prev(i);
      // get initial velocity of the curve after this angle
      gsMatrix<T> velNext = trimSurface.TangentCoefs_next(i);

      // compute angle
      angles.push_back(math::atan2(velNext(1), velNext(0)) - math::atan2(velPrev(1), velPrev(0)));
      if(angles[i] <= -EIGEN_PI) angles[i] += 2 * EIGEN_PI;
      if(angles[i] > EIGEN_PI) angles[i] -= 2 * EIGEN_PI;
      GISMO_ASSERT(angles[i] != EIGEN_PI, "Turning angle of pi is invalid");

      // replace convex angles by straight lines
      if(angles[i] < 0) angles[i] = 0;

      // get distance to the next angle
      lengths.push_back((trimSurface.vertexCoord(0, (i + 1) % n) - trimSurface.vertexCoord(0, i)).norm());
    }
  }
  gsCurveLoop<T>::approximatingPolygon(angles, lengths, 0.1, corners);

}

template<class T>
void gsVolumeSegment<T>::attemptTraceCurve(gsTrimSurface<T> &trimSurface, const gsMatrix<T> cornerPoints,
                                           gsTemplate<T> &templ, const gsMatrix<T> &firstTrialPoint, const gsMatrix<T> &step,
                                           int maxTrials, bool harmonic, int n_fit_points, T traceTolerance, gsMatrix<T> &result)
{
  std::pair<gsFunction<T>*, gsFunction<T>*> components;
  gsMVInterpolation<T> * mvi = 0;
  if(harmonic)
  {
    // use harmonic functions

      GISMO_ERROR("harmonic mapping used here.");
//    components = trimSurface.domain().mapto(templ);
  }
  else
  {
    // use mean value interpolation
    mvi = new gsMVInterpolation<T>(&trimSurface, cornerPoints);
    components.first = new gsMVInterpolationComponent<T>(mvi, 0);
    components.second = new gsMVInterpolationComponent<T>(mvi, 1);
  }
  int trials = 0;
  bool done;
  gsMatrix<T> start;
  do
  {
    // TODO: try a new value for start
    start = firstTrialPoint + step * trials;
    gsTraceCurve<T>(components, start, templ.skeleton(0), result, n_fit_points, traceTolerance);
    trials++;
    done = true; // assume done and adjust by looking for infinite values
    for(index_t i = 0; i < result.rows(); i++)
    {
      done = done && !(math::isnan(result(i, 0)) || math::isnan(result(i, 1)));
    }
  } while(trials < maxTrials && !done);
  GISMO_ASSERT(done, "Curve tracing failed to converge for all starting locations");
  if(!harmonic)
  {
    delete components.first;
    delete components.second;
    delete mvi;
  }
}

template<class T>
void gsVolumeSegment<T>::chooseCuttingCurve(gsTrimSurface<T> &trimSurface,index_t idxStart, index_t idxEnd,
                                            T w_reg, T w_app, int n_fit_points, T traceTolerance,
                                            bool harmonic, bool regular, T curveFraction, int interpDegree,
                                            gsBSpline<T> **imgCurve,
                                            gsBSpline<T> &cuttingCurve, gsMatrix<T> &fitPoints)
{
  const gsCurveLoop<T> & curveLoop = trimSurface.domain().outer();
  const index_t n = (index_t)curveLoop.size();

  // endpoints form the exact constraints
  gsMatrix<T> boundaryTimes(1, 2);
  boundaryTimes << 0, 1;
  gsMatrix<T> boundaryValues(2, 2);
  boundaryValues.col(0) = curveLoop.curve(idxStart).coefs().row(0);
  boundaryValues.col(1) = curveLoop.curve(idxEnd).coefs().row(0);

  // compute bisectors at the endpoints, these form additional constraints
  gsMatrix<T> corJacobianStart = trimSurface.derivatives(idxStart);
  gsMatrix<T> corJacobianEnd = trimSurface.derivatives(idxEnd);
  gsMatrix<T> velStartPrev = trimSurface.UnitTangentCoefs_prev(idxStart, corJacobianStart);
  gsMatrix<T> velStartNext = trimSurface.UnitTangentCoefs_next(idxStart, corJacobianStart);
  gsMatrix<T> velEndPrev = trimSurface.UnitTangentCoefs_prev(idxEnd, corJacobianEnd);
  gsMatrix<T> velEndNext = trimSurface.UnitTangentCoefs_next(idxEnd, corJacobianEnd);
  gsMatrix<T> normalValues(2, 2);
  T nvsn = velStartNext.squaredNorm();
  T nvsp = velStartPrev.squaredNorm();
  T nven = velEndNext.squaredNorm();
  T nvep = velEndPrev.squaredNorm();
  velStartNext *= nvsp;
  velStartPrev *= nvsn;
  velEndNext *= nvep;
  velEndPrev *= nven;
  // (bisectors will be orthogonal to these differences)
  normalValues.col(0) = velStartNext - velStartPrev;
  normalValues.col(1) = velEndNext - velEndPrev;

  gsMatrix<T> cornerPoints;
  gsVolumeSegment<T>::chooseConvexPolygon(trimSurface, regular, cornerPoints);

  // choose an interval that we will search for a good starting point. we do this
  // by solving where the boundary intersects the bisector of the angle at the
  // starting vertex. This is a union of intervals along the bisector, and one
  // of these intervals has the starting vertex as an endpoint. that's the
  // interval we choose.
  gsMatrix<T> m = (velStartNext - velStartPrev).transpose();
  gsMatrix<T> v(1, 2);
  v << -m(0, 1), m(0, 0);
  typename gsBSpline<T>::uPtr boundarySingleCurve = memory::convert_ptr<gsBSpline<T> >(trimSurface.domain().loop(0).singleCurve());
  GISMO_ASSERT(boundarySingleCurve != NULL, "Only implemented for B-splines");
  gsMatrix<T> boundaryCoefs = boundarySingleCurve->coefs();
  gsMatrix<T> coefs1d(boundaryCoefs.rows(), 1);
  gsMatrix<T> p = boundaryValues.col(0).transpose();
  for(index_t i = 0; i < coefs1d.rows(); i++)
  {
    coefs1d(i, 0) = ((boundaryCoefs.row(i) - p) * m.transpose()).value(); // row vector times column vector
  }
  gsBSpline<T> orthog(boundarySingleCurve->basis(), give(coefs1d));

  // find all the intersections between a straight line and the boundary.
  // slv.allRoots will return the roots as parameters of the boundary curve
  std::vector<T> intersectionTimes;
//  gsBSplineSolver<T> slv;
//  slv.allRoots(orthog, intersectionTimes);

  // FIXME this is a temporary fix.
  // Both gsBSplineSolver and findHyperPlaneIntersections are known to be buggy.
  // The second can fail in some cases of tangential intersections, but it is generally more reliable.
  // Also part of the following could be simplified as findHyperPlaneIntersections gives more information.
  // See also the todo at line 402.
  std::vector<Root<T> > roots;
  gsVector<T> normal;
  normal.setConstant(1,1);
  unsigned numRoots = findHyperPlaneIntersections<T>(orthog,normal, 0, 100*std::numeric_limits<T>::epsilon(), roots);
  for(unsigned i=0; i<numRoots;++i)
  {
    intersectionTimes.push_back(roots[i].parameter);
  }

  // if the first intersection is at the start of the boundary curve then it
  // will be repeated as the last vertex, since the curve is closed. we need
  // to remove it.
  gsMatrix<T> boundaryCurveRange = boundarySingleCurve->parameterRange();
  if(intersectionTimes[0] <= boundaryCurveRange(0, 0) &&
     intersectionTimes[intersectionTimes.size() - 1] >= boundaryCurveRange(0, 1))
  {
    intersectionTimes.resize(intersectionTimes.size() - 1);
  }
  // convert each intersection from a parameter of the B-spline to a
  // parameter of the line
  gsMatrix<T> t(1, 1), x;
  for(size_t i = 0; i < intersectionTimes.size(); i++)
  {
    t(0,0) = intersectionTimes[i];
    boundarySingleCurve->eval_into(t, x);
    intersectionTimes[i] = (v * (x - p.transpose())).value() / (v * v.transpose()).value();
  }

  // Find the first positive root with an odd index.
  // TODO: I think we could find the interval without sorting the points and then
  // looping through them all... just have to be aware that numerical error could
  // cause the root corresponding to the start vertex to be nonzero.
  T exitTime = T(-1);
  std::sort(intersectionTimes.begin(), intersectionTimes.end());
  for(size_t i = 1; i < intersectionTimes.size(); i+=2)
  {
    if(intersectionTimes[i] > 0)
    {
      exitTime = intersectionTimes[i];
      break;
    }
  }
  GISMO_ASSERT(exitTime > 0, "Failed to find intersection of bisector with outer loop");

  // create a template with the curve to be traced.
  std::vector<T> domainSizes = curveLoop.domainSizes();
  gsTemplate<T> templ(domainSizes, cornerPoints, idxStart, idxEnd, curveFraction);
  if(imgCurve != NULL) *imgCurve = (gsBSpline<T>*)templ.skeleton(0)->clone().release();
  std::vector<T> intervalWidths;
  for(index_t i = 0; i < n; i++)
  {
    gsMatrix<T> range = curveLoop.curve(i).parameterRange();
    GISMO_ASSERT(range.rows() == 1 && range.cols() == 2, "Expected 1x2 matrix for parameter range of curve");
    intervalWidths.push_back(range(0, 1) - range(0, 0));
  }
  gsMatrix<T> traceTimes =  gsPointGrid<T>((T(1) - curveFraction) / 2, (T(1) + curveFraction) / 2, n_fit_points);
  // evaluate the points and check that they are inside the convex polygon
  // TODO: this problem could be prevented by forcing all the control points
  // to stay inside the convex polygon.
  GISMO_ASSERT(cornerPoints.rows() == n, "Unexpected number of corner points");
  for(index_t j = 0; j < n; j++)
  {
    gsMatrix<T> s = (cornerPoints.row((j + 1) % n) - cornerPoints.row(j));
    gsMatrix<T> srot(2, 1);
    srot << -s(0, 1), s(0, 0);
    for(int i = 0; i < traceTimes.cols(); i++)
    {
      gsMatrix<T> d = templ.skeleton(0)->eval(traceTimes.col(i)).transpose();
      gsMatrix<T> crossProd = (d - cornerPoints.row(j)) * srot;
      GISMO_ASSERT(crossProd.size() == 1, "Unexpected matrix product");
      GISMO_ASSERT(crossProd(0, 0) >= 0, "Curve to trace is not inside polygon - treatment of this case has not been implemented");
    }
  }
  // trace curve
  // (this bisection template only has one skeleton curve in it: templ.skeleton(0).)
  int numTrials = 10;
  gsMatrix<T> trialStep = v.transpose() * exitTime / (numTrials + 1);
  gsMatrix<T> tracedPts;
  attemptTraceCurve(trimSurface, cornerPoints, templ, boundaryValues.col(0) + trialStep, trialStep,
                    numTrials, harmonic, n_fit_points, traceTolerance, tracedPts);

  gsMatrix<T> support = trimSurface.basis().support();
  // what if some of the points found in attemptTraceCurve didn't converge?
  // clean them out before we fit the curve.
  size_t idxOut = 0;
  GISMO_ASSERT(tracedPts.rows() == 2, "Unexpected dimension for traced points");
  for(index_t idxIn = 0; idxIn < tracedPts.cols(); idxIn++)
  {
    bool ok = true;
    for(index_t dim = 0; dim < 2; dim++)
    {
      if(tracedPts(dim, idxIn) < support(dim, 0)||
         tracedPts(dim, idxIn) > support(dim, 1)) ok = false;
    }
    if(ok)
    {
      tracedPts.col(idxOut) = tracedPts.col(idxIn);
      traceTimes(0, idxOut) = traceTimes(0, idxIn);
      idxOut++;
    }
  }
  gsMatrix<T> preimageApp = traceTimes.block(0, 0, 1, idxOut);
  fitPoints = tracedPts.block(0, 0, 2, idxOut);

  // create a spline
  gsKnotVector<T> kv(0., 1., 0, interpDegree + 1);
  gsMatrix<T> pointResiduals, normalResiduals;
  cuttingCurve = gsInterpolate(kv, boundaryTimes, boundaryValues, boundaryTimes, normalValues,
                        preimageApp, fitPoints, w_reg, w_app, pointResiduals, normalResiduals);
  // check if numerical error has moved the endpoints the support of
  // the master surface, and adjust if necessary
  const T endPointTolerance = 4 * 2.220446049250313e-16;
  GISMO_UNUSED(endPointTolerance);

  gsMatrix<T> &coefs = cuttingCurve.coefs();
  for(index_t i = 0; i < coefs.rows(); i += coefs.rows() - 1) // only check the endpoints
  {
    for(index_t dim = 0; dim < 2; dim++)
    {
      assert(coefs(i, dim) >= support(dim, 0) - endPointTolerance);
      coefs(i, dim) = math::max(coefs(i, dim), support(dim, 0));
      assert(coefs(i, dim) <= support(dim, 1) + endPointTolerance);
      coefs(i, dim) = math::min(coefs(i, dim), support(dim, 1));
    }
  }
  gsMatrix<T> dotProd = (coefs.row(1) - coefs.row(0)) * normalValues.col(0);
  GISMO_ASSERT(math::abs(dotProd(0, 0)) < 0.0001, "Normal condition not satisfied");
  dotProd = (coefs.row(coefs.rows() - 1) - coefs.row(coefs.rows() - 2)) * normalValues.col(1);
  GISMO_ASSERT(math::abs(dotProd(0, 0)) < 0.0001, "Normal condition not satisfied");

//  for(index_t i = 0; i < coefs.rows(); i++)
//  {
//    for(index_t dim = 0; dim < 2; dim++)
//    {
//      assert(coefs(i, dim) >= support(dim, 0) - endPointTolerance);
//      assert(coefs(i, dim) <= support(dim, 1) + endPointTolerance);
//    }
//  }

  gsMatrix<T> innerBisectors(2, 2);
  // assert that the direction of the spline at the endpoints points inwards
  innerBisectors.col(0) << -normalValues(1, 0), normalValues(0, 0);
  innerBisectors.col(1) << -normalValues(1, 1), normalValues(0, 1);
  assert(((coefs.row(1) - coefs.row(0)) * innerBisectors.col(0))(0, 0) > 0);
  assert(((coefs.row(coefs.rows() - 2) - coefs.row(coefs.rows() - 1)) * innerBisectors.col(1))(0, 0) > 0);
}


template<class T>
void gsVolumeSegment<T>::splitFaces(gsSolid<T> &sl, const gsCuttingLoop<T> &cuttingLoop,
                                    T w_reg, T w_app, int n_fit_points, T traceTolerance,
                                    bool harmonic, T curveFraction)
{
  std::vector<unsigned> path = cuttingLoop.computePath();
  std::vector<unsigned>::size_type numVerts = path.size();
  gsMatrix<T> fitPoints;
  for(std::vector<unsigned>::size_type i = 0; i < numVerts; i++)
  {
    gsSolidHeVertex<T> *vert1 = sl.vertex[path[i]];
    gsSolidHeVertex<T> *vert2 = sl.vertex[path[(i + 1) % numVerts]];
    
    if(vert1->hasHalfEdge(vert2))
    {
      // edge already exists, no need to create
      continue;
    }
    
    // edge does not exist, need to create it. first find the face containing both vertices.
    std::vector<gsSolidHalfFace<T> *> faces = vert1->getHalfFaces();
    typename std::vector<gsSolidHalfFace<T> *>::size_type numFaces = faces.size();
    gsSolidHalfFace<T> * containingFace = NULL;
    for(typename std::vector<gsSolidHalfFace<T> *>::size_type j = 0; j < numFaces; j++)
    {
      if(faces[j]->containsVertex(vert2))
      {
        containingFace = faces[j];
      }
    }
    GISMO_ASSERT(containingFace != NULL, "must find a face containing both vertices.");
    
    // find the indices of these vertices within the face.
    int idxStart = containingFace->indexOfVertex(vert1);
    int idxEnd = containingFace->indexOfVertex(vert2);
    gsTrimSurface<T> *trimSurface = containingFace->surf;
    gsBSpline<T> *auxE = new gsBSpline<T>;
    gsVolumeSegment<T>::chooseCuttingCurve(*trimSurface, idxStart, idxEnd, w_reg, w_app,
                                           n_fit_points, traceTolerance, harmonic, false, curveFraction,
                                           5, NULL, *auxE, fitPoints);

    // split this face along the spline.
    sl.splitFace(containingFace, vert1, vert2, auxE);
  }

  sl.checkStructure();
}


template <class T>
gsTrimSurface<T>* gsVolumeSegment<T>::cuttingSurface(  gsCuttingLoop<T> const & cloop, const gsCurveLoop<T> & curveLoop, gsKnotVector<T> const & kv1,
                  gsKnotVector<T> const & kv2,gsInterpOption<T> & intpOpt, bool force_normal) const
{
    const int nInnerPoints = intpOpt.ndInPoint;
    const T wPoint = intpOpt.wdPoint;
    const T wNormal = intpOpt.wdNormal;
    const T wEnergy = intpOpt.wReg;
    const int nPoints=2+nInnerPoints;

    GISMO_ASSERT( static_cast<size_t>(nPoints) >= kv1.size() - kv1.degree() - 1 && 
                  static_cast<size_t>(nPoints) >= kv2.size() - kv2.degree() - 1, 
                  "Not enough points to use this knot vector");

    //- get coordinates of corners of the cutting loop
    gsMatrix<T> corners3D = cloop.getVertices();
    //- get coordinates of corners in the parameter domain
    gsMatrix<T> corners2D = curveLoop.sample(2,1);

    // eval splines at parameter values
    gsMatrix<T> bInner2D = curveLoop.sample(nPoints,2);
    gsMatrix<T> bInner3D = cloop.sample(nPoints,2);

    // eval normals
    gsMatrix<T> normalv = cloop.sampleNormal(nInnerPoints+2);
    gsMatrix<T> dummy2D(2, 0), dummy3D(3, 0);

    typename gsTensorBSpline<2,T>::Ptr master = gsInterpolateSurface(corners2D, corners3D,
                       bInner2D, bInner3D,
                       dummy2D, dummy3D,
                       bInner2D, normalv,
                       wPoint, T(0), wNormal, wEnergy,
                       kv1, kv2, force_normal);

    gsCurveLoop<T> * newCurveLoop = new gsCurveLoop<T>(curveLoop);
    gsPlanarDomain<T> * domain= new gsPlanarDomain<T>(newCurveLoop);

    return new gsTrimSurface<T>(master,domain);
}


} // namespace 

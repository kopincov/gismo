/** @file gsMVInterpolation.hpp

    @brief Provides implementation of gsMVInterpolation class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Pauley
*/

#pragma once

#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsTraceCurve.hpp>

namespace gismo
{


/// compute an integral needed for mean value interpolation.
/// \a B1 and \a coeffs are the basis and coefficients of the spline,
/// \a pts is a matrix whose columns are the points you want to interpolate
/// at, and \a val1 and \a val2 are the values at the endpoints of the
/// spline of the function you want to interpolate.
template<class T>
typename gsMatrix<T>::uPtr MVCint( const gsBasis<T>& B1, 
                                   const gsMatrix<T>& coefs, 
                                   const gsMatrix<T>& pts, 
                                   const gsVector<T>& val1, 
                                   const gsVector<T>& val2, int nGauss )
{
    index_t numPts = pts.cols();
    index_t targetDim = val1.size();

    assert(coefs.cols() == 2);
    assert(pts.rows() == 2);
    assert(targetDim == val2.size());

    std::vector<T> breaks = B1.domain()->breaks();
    gsMatrix<T> ngrid;
    gsVector<T> wgrid;
    // Quadrature points
    gsGaussRule<T>(nGauss).mapToAll(breaks, ngrid, wgrid);
  
    gsMatrix<T> range = B1.support();
    T mint = range(0, 0);
    T maxt = range(0, 1);
    
    gsMatrix<T>  ev  = B1.eval(ngrid); // Evaluate over the grid
    gsMatrix<T>   evDeriv  = B1.deriv(ngrid); // Evaluate over the grid
  
    gsMatrix<T> * integral = new gsMatrix<T>(targetDim, numPts);
    gsMatrix<T> integrand(targetDim, numPts);
    integral->setZero();
    gsVector<T,2> c, cDeriv;
  
    for (index_t k=0; k!= ngrid.cols(); ++k)
    {
        c.setZero();
        cDeriv.setZero();
        for (index_t i=0; i != ev.rows(); ++i)
        {
            c += ev(i, k) * coefs.template block<1,2>(i,0);
            cDeriv += evDeriv(i, k) * coefs.template block<1,2>(i,0);
        }
    
        gsVector<T> interpVal = ((maxt - ngrid(k)) * val1 + (ngrid(k) - mint) * val2) / (maxt - mint);
        for(index_t thisCol = 0; thisCol < numPts; thisCol++)
        {
            T wNumerator = (pts(1, thisCol) - c(1)) * cDeriv(0) - (pts(0, thisCol) - c(0)) * cDeriv(1);
            T wDenominator = math::pow((pts.col(thisCol) - c).squaredNorm(), 1.5);
            integrand.col(thisCol) = interpVal * (wNumerator / wDenominator);
        }
        *integral += wgrid(k) * integrand;
    }

    return typename gsMatrix<T>::uPtr(integral);
}


/// Compute the derivative of MVCint. Returns a (n x (2 m)) sized
/// matrix (where n == val1.size() == val2.size(), and m = pts.cols())
/// with the Jacobians stored side-by-side (to agree with the
/// convention of other classes such as gsGeometry).
template<class T>
typename gsMatrix<T>::uPtr MVCintGrad( const gsBasis<T>& B1, 
                                       const gsMatrix<T>& coefs, 
                                       const gsMatrix<T>& pts, 
                                       const gsVector<T>& val1, 
                                       const gsVector<T>& val2, int nGauss )
{
    index_t targetDim = val1.size();
    index_t numPts = pts.cols();

    assert(coefs.cols() == 2);
    assert(pts.rows() == 2);
    assert(targetDim == val2.size());

    std::vector<T> breaks = B1.domain()->breaks();
    gsMatrix<T> ngrid;
    gsVector<T> wgrid;
    // Quadrature points
    gsGaussRule<T>(nGauss).mapToAll(breaks, ngrid, wgrid);

    gsMatrix<T> range = B1.support();
    T mint = range(0, 0);
    T maxt = range(0, 1);
    
    gsMatrix<T>   ev  = B1.eval(ngrid); // Evaluate over the grid
    gsMatrix<T>   evDeriv  = B1.deriv(ngrid); // Evaluate over the grid
  
    gsMatrix<T> * integral = new gsMatrix<T>(targetDim, 2 * numPts);
    gsMatrix<T> gradAtPt(targetDim, 2 * numPts);
    integral->setZero();
    gsVector<T,2> c, cDeriv, rotDeriv;

    for (index_t k=0; k!= ngrid.cols(); ++k)
    {
        c.setZero();
        cDeriv.setZero();
        for (index_t i=0; i != ev.rows(); ++i)
        {
            c += ev(i, k) * coefs.template block<1,2>(i,0);
            cDeriv += evDeriv(i, k) * coefs.template block<1,2>(i,0);
        }
    
        gsVector<T> interpVal = ((maxt - ngrid(k)) * val1 + (ngrid(k) - mint) * val2) / (maxt - mint);
        rotDeriv << -cDeriv(1), cDeriv(0);
        for(index_t thisCol = 0; thisCol < numPts; thisCol++)
        {
            gsVector<T,2> disp = c - pts.col(thisCol);
            T sqNormDisp = disp.squaredNorm();
            gradAtPt.block(0, 2 * thisCol, targetDim, 2) =
                interpVal * ((sqNormDisp * rotDeriv) - (3 * disp.dot(rotDeriv) * disp)).transpose() /
                math::pow(sqNormDisp, 2.5);
        }
    
        *integral += wgrid(k) * gradAtPt;
    }

    return typename gsMatrix<T>::uPtr(integral);
}


template<class T>
void gsMVInterpolation<T>::init(int nGauss2)
{
    this->cached = false;
    this->nGauss = nGauss2;

    assert(this->trimSurface->domain().outer().size() == (int)interpolationPoints.rows());
    assert(this->interpolationPoints.cols() == 2);
}

template<class T>
void gsMVInterpolation<T>::updateCalculations(const gsMatrix<T> *u) const
{
    if(cached && (u->rows() == cachedEvalPts.rows() && u->cols() == cachedEvalPts.cols() && *u == cachedEvalPts))
    {
        return;
    }
    const gsCurveLoop<T> & curveLoop = this->trimSurface->domain().outer();
    size_t numCurves = curveLoop.size();
    size_t numEvalPts = u->cols();

    assert(u->rows() == 2);
    assert(numCurves == (size_t)interpolationPoints.rows());

    gsMatrix<T> weightIntegral(1, numEvalPts);
    weightIntegral.setZero();
    gsMatrix<T> dw(1, 2 * numEvalPts);
    dw.setZero();
    gsMatrix<T> fw(2, numEvalPts);
    fw.setZero();

    gsVector<T> one(1);
    one << 1;

    // 2x2 Jacobians are to be stored in blocks: [ J1  J2  J3 ... Jn ]
    gsMatrix<T> fdw(2, 2 * numEvalPts);
    fdw.setZero();

    gsVector<T> ip1, ip2;

    for(size_t i = 0; i < numCurves; i++)
    {
        const gsGeometry<T> & thisCurve = curveLoop.curve(i);

        ip1 = interpolationPoints.row(i);
        ip2 = interpolationPoints.row((i + 1) % numCurves);

        fw += *MVCint(
            thisCurve.basis(), thisCurve.coefs(),
            *u, ip1, ip2, this->nGauss);
        fdw += *MVCintGrad(
            thisCurve.basis(), thisCurve.coefs(),
            *u, ip1, ip2, this->nGauss);
        weightIntegral += *MVCint(thisCurve.basis(), thisCurve.coefs(), *u, one, one, this->nGauss);
        dw += *MVCintGrad(thisCurve.basis(), thisCurve.coefs(), *u, one, one, this->nGauss);
    }

    cached = true;
    cachedEvalPts = *u;
    cachedValue.resize(2, numEvalPts);
    cachedDeriv.resize(2, 2 * numEvalPts);
    for(size_t thisPt = 0; thisPt < numEvalPts; thisPt++)
    {
        cachedValue.col(thisPt) = fw.col(thisPt) / weightIntegral(0, thisPt);
        cachedDeriv.template block<2,2>(0, thisPt * 2) =
            (fdw.template block<2,2>(0, thisPt * 2) - 
             cachedValue.template block<2,1>(0, thisPt) * dw.template block<1,2>(0, thisPt * 2))
            / weightIntegral(0, thisPt);
    }

}

template<class T>
gsMatrix<T> gsMVInterpolation<T>::support() const
{
    return this->trimSurface->domain().outer().getBoundingBox();
}

template<class T>
void gsMVInterpolation<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    updateCalculations(&u);
    result = cachedValue;
}

template<class T>
void gsMVInterpolation<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    updateCalculations(&u);
    result = cachedDeriv;
}

template<class T>
void gsMVInterpolation<T>::eval_component_into(const gsMatrix<T>& u, 
                                               const index_t comp, 
                                               gsMatrix<T>& result) const
{
    updateCalculations(&u);
    result = cachedValue.row(comp);
}

template<class T>
void gsMVInterpolation<T>::deriv_component_into(const gsMatrix<T>& u, 
                                                const index_t comp, 
                                                gsMatrix<T>& result) const
{
    updateCalculations(&u);
    result = cachedDeriv.row(comp);
}

template<class T>
gsMatrix<T> gsMVInterpolationComponent<T>::support() const
{
    return parent->support();
}

template<class T>
void gsMVInterpolationComponent<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    return parent->eval_component_into(u, idx, result);
}

template<class T>
void gsMVInterpolationComponent<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    return parent->deriv_component_into(u, idx, result);
}

}// namespace gismo

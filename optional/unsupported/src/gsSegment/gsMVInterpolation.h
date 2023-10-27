/** @file gsMVInterpolation.h

    @brief Provides declaration of gsMVInterpolation class for computing mean
    value interpolation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Pauley
*/

#pragma once

#include <gsCore/gsFunction.h>

namespace gismo
{

template<class T> class gsMVInterpolationComponent;

/// @brief Class for evaluating the mean value interpolation mapping the interior of
/// a curve loop to the interior of a convex polygon. 
///
/// Since a lot of the
/// work is common between computing position and the Jacobian, this class
/// tries to save some effort by remembering the last set of points it was
/// asked about and the value of the function and derivative there.
template<class T>
class gsMVInterpolation : public gsFunction<T>
{

public:
    /// Shared pointer for gsMVInterpolation
    typedef memory::shared_ptr< gsMVInterpolation > Ptr;

    /// Unique pointer for gsMVInterpolation
    typedef memory::unique_ptr< gsMVInterpolation > uPtr;

    /// Construct a mean value interpolation from a given curve loop and the
    /// values of the image at the end points of the curves.
    /// \param _trimSurface The trimmed surface to use as the domain
    /// \param _interpolationPoints the corners of the image
    /// \param nGauss number of Gauss points to use for integration
    gsMVInterpolation(gsTrimSurface<T> * _trimSurface, 
                      const gsMatrix<T> & _interpolationPoints, 
                      int nGauss = 32)
    {
        this->trimSurface = _trimSurface;
        this->interpolationPoints = _interpolationPoints; // copy
        init(nGauss);
    }

    /// Look at a point and update calculations.
    /// \param u Point to look at.
    void updateCalculations(const gsMatrix<T> *u) const;

    virtual short_t domainDim () const {return 2;}
    virtual short_t targetDim () const {return 2;}
    
    // overloads of virtual functions
    GISMO_CLONE_FUNCTION(gsMVInterpolation)
    virtual gsMatrix<T> support() const;
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    virtual void eval_component_into(const gsMatrix<T>& u, 
                                     const index_t comp, 
                                     gsMatrix<T>& result) const;

    virtual void deriv_component_into(const gsMatrix<T>& u, 
                                      const index_t comp, 
                                      gsMatrix<T>& result) const;

        private:

    /// internal initializer, used by the constructors
        void init(int nGauss);

    /// choose interpolation points based on cutting loop
    void automaticInterpolationPoints()
    { GISMO_NO_IMPLEMENTATION };

    gsTrimSurface<T> *trimSurface;
    gsMatrix<T> interpolationPoints;

    mutable bool cached;
    /// Stores the matrix of points this function was last evaluated at.
    mutable gsMatrix<T> cachedEvalPts;

    /// Stores the matrix whose columns are the values of the function evaluated
    /// at the points given by the columns of cachedEvalPts.
    mutable gsMatrix<T> cachedValue;

    /// Stores the 2x2 Jacobians in blocks: [ J1  J2  J3 ... Jn ].
    mutable gsMatrix<T> cachedDeriv;

    /// Number of points to use for Gaussian quadrature
    int nGauss;
};

///Class representing one component of the mean value interpolation
// to do: replace
template<class T>
class gsMVInterpolationComponent : public gsFunction<T>
{
public:
    /// Shared pointer for gsMVInterpolationComponent
    typedef memory::shared_ptr< gsMVInterpolationComponent > Ptr;

    /// Unique pointer for gsMVInterpolationComponent
    typedef memory::unique_ptr< gsMVInterpolationComponent > uPtr;

    /// constructor
    gsMVInterpolationComponent(gsMVInterpolation<T> * _parent, index_t _idx)
    {
        this->parent = _parent;
        this->idx    = _idx;
    }

    // overloads of virtual functions
    GISMO_CLONE_FUNCTION(gsMVInterpolationComponent)
    virtual gsMatrix<T> support() const;
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    virtual short_t domainDim () const {return 2;}
    virtual short_t targetDim () const {return 1;}

private:

    gsMVInterpolation<T> * parent;
    index_t idx;

};


}


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMVInterpolation.hpp)
#endif


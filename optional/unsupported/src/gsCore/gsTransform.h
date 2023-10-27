/** @file gsTransform.h

    @brief Represents an isogeometric geometry map

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo {

/**
   @brief this object encaptulates the computation needed to transform
   values/derivatives from a parametric domain to a physical domain

   Derived from this class are: 
   - gsTransformId
   - gsTransformHgrad, 
   - gsTransformHdiv
   - gsTransformHcurl

 */
template <typename T>
class gsTransform 
{
protected:
    typedef typename gsMatrix<T>::constColumn  constColumn ;
    typedef typename gsMatrix<T>::constColumns constColumns;

public:
    
    gsTransform() : points(NULL) //, map(NULL)
    { }

protected:

    // Pointer to a set of points in parameter domain of \a map
    const gsMatrix<T> * points;

    // Pointer to a geometry map function (probably not needed)
    //const gsFunction<T> * map;
    
    // Contains computed values of the \a map
    gsFuncData<T> m_data;
    
    // Extra computed quantities related to the \a map
    unsigned m_extraFlags;

    gsMatrix<T> m_measures;
    gsMatrix<T> m_gradTransforms;
    gsMatrix<T> m_normals;
    
    unsigned m_parDim;

public:

    const gsFuncData<T> & data() const { return m_data;}

    const gsMatrix<T> & values() const
    { 
        return m_data.values[0];
    }

    const gsMatrix<T> & params() const
    {
        GISMO_ASSERT( points!=NULL, "No parameter points given to gsTransform");
        return *points;
    }

    const gsMatrix<T> & measures() const 
    { 
        GISMO_ASSERT(m_extraFlags| NEED_MEASURE, 
                     "measures are not computed unless the NEED_MEASURE flag is set."); 
        return m_measures;
    }

    const constColumn measure(const index_t k) const 
    { 
        GISMO_ASSERT(m_extraFlags| NEED_MEASURE, 
                     "measures are not computed unless the NEED_MEASURE flag is set."); 
        return m_gradTransforms.col(k);
    }

    const gsMatrix<T> & gradTransforms() const 
    { 
        GISMO_ASSERT(m_extraFlags| NEED_GRAD_TRANSFORM, 
                     "gradient transformation matrice are not computed unless the NEED_MEASURE flag is set."); 
        return m_gradTransforms;
    }
    
    const constColumns gradTransform(const index_t k) const 
    { 
        GISMO_ASSERT(m_extraFlags| NEED_GRAD_TRANSFORM, 
                     "gradient transformation matrice are not computed unless the NEED_MEASURE flag is set."); 
        return m_gradTransforms.middleCols( k*m_parDim, m_parDim);
    }

    const gsMatrix<T> & normals() const 
    { 
        GISMO_ASSERT(m_extraFlags| NEED_NORMAL, 
                     "normals are not computed unless the NEED_NORMAL flag is set."); 
        return m_normals;
    }

    const constColumn normal(const index_t k) const 
    { 
        GISMO_ASSERT(m_extraFlags| NEED_GRAD_TRANSFORM, 
                     "normals are not computed unless the NEED_NORMAL flag is set."); 
        return m_normals.col(k);
    }

public:

    ///\brief Computes quantities implied by flags on the input
    ///parametric \a a_points through the transformation defined by \a
    ///_map
    void compute(const gsMatrix<T> & _points, const gsFunction<T> & _map) //, unsigned flags
    {
        points   = &points;
        // map   = &map;
        m_parDim = _map.domainDim();

        _map.compute(_points, m_data);

        // Based on flags:
        // Compute measure
        //m_measures = 

        // Compute Jacobian inverses
        // ..

        // Compute normals
        // ..
    }

    
    /// Assumes that \dataIn are values at the same points where
    /// compute(.,.) was called on
    virtual void transform(const gsFuncData<T> & dataIn,
                           gsFuncData<T> & result) const = 0;
    /*
    {
        if ( ! dataIn.isParametric() )
            return;

        // Apply transformation on the quantities (values, gradients,
        // 2nd ders..) in dataIn and store the transformed values in dataOut

    }
    */

    /// Assumes that \dataInOut are values at the same points where
    /// compute(.,.) was called on
    virtual void transformInPlace(gsFuncData<T> & dataInOut)  const
    {
        //if ( ! dataInOut.isParametric() )
        //    return;

        gsFuncData<T> tmp;
        transform(dataInOut,tmp);
        tmp.swap(dataInOut);

        // Apply transformation on the quantities (values, gradients,
        // 2nd ders..) in-place        
    }

};


template <typename T>
class gsTransformId : public gsTransform<T>
{
public:
    //problem: values() will not work ? unless we re-impl. compute(.,.)

    void transform(const gsFuncData<T> & dataIn,
                   gsFuncData<T> & result) const
    {
        // Pass on the input to output
        result = dataIn;
    }

    void transformInPlace(gsFuncData<T> & dataInOut)  const
    {
        // Do nothing
    }
};


template <typename T>
class gsTransformHgrad : public gsTransform<T>
{
public:

    void transform(const gsFuncData<T> & dataIn,
                   gsFuncData<T> & result) const
    {
        const unsigned flags = dataIn.flags;
        
        if (flags & NEED_VALUE)
        {
            result.values.resize(1);
            result.values[0] = dataIn.values[0];
        }
/*
        if ( flags & NEED_DERIV ) // note: NEED_GRAD
        {
            result.values.resize(2);
            GISMO_ASSERT(m_flags & NEED_GRAD_TRANSFORM, "J^-1 not computed");
            GISMO_ASSERT(allGrads.rows() % ParDim == 0, "Invalid size of gradient matrix");

            const index_t numGrads = allGrads.rows() / ParDim;
            const gsAsConstMatrix<T,ParDim> grads_k(allGrads.col(k).data(), ParDim, numGrads);
            trfGradsK.noalias() = m_jacInvs.template block<GeoDim,ParDim>(0, k*ParDim) * grads_k;


            result.values[1] = gradTransforms() * dataIn.values[1];

            deriv_into(points,result.values[1]);
            if( flags & NEED_DIV )
                convertValue<T>::derivToDiv(result.values[1], result.divs, map->info());
            if( flags & NEED_CURL )
                convertValue<T>::derivToCurl(result.values[1], result.curls, map->info());
        }
        else
        {
            if (flags & NEED_DIV)
                div_into(points, result.divs);
            if (flags & NEED_CURL)
                curl_into(points, result.curls);
        }
        
        if (flags & NEED_DERIV2)
        {
            deriv2_into(points, result.values[2]);
            if (flags & NEED_LAPLACIAN)
                convertValue<T>::deriv2ToLaplacian(result.values[2], result.laplacian, map->info());
        }
        else
        {
            if (flags & NEED_LAPLACIAN)
                laplacian_into(points, result.laplacian);
        }

//*/
    }

    void transformInPlace(gsFuncData<T> & dataInOut)  const
    {
        //const unsigned flags = dataInOut.flags;

    }

};


}//namespace gismo


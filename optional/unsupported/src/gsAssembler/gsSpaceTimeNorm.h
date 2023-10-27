/** @file gsSpaceTimeNorm.h

    @brief Computes the Space-Time norm.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Moore
*/

#include<gsAssembler/gsNorm.h>
#include <gsCore/gsBoundary.h>
#include <gsCore/gsBoxTopology.h>


#pragma once

namespace gismo
{

/** @brief The gsSpaceTimeNorm class provides the functionality
 * to calculate the Space-Time norm between a field and a function.
 *
 * \ingroup Assembler
*/
template <class T>
class gsSpaceTimeNorm : public gsNorm<T>
{
    friend  class gsNorm<T>;
    typedef typename gsMatrix<T>::RowsBlockXpr Rows;

public:

    gsSpaceTimeNorm(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             bool _f2param = false) 
    : gsNorm<T>(_field1,_func2), dfunc2(NULL), f2param(_f2param)
    { 
        m_theta = 1;
    }

    gsSpaceTimeNorm(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             const gsFunction<T> & _dfunc2,
             bool _f2param = false)
    : gsNorm<T>(_field1,_func2), dfunc2(&_dfunc2), f2param(_f2param)
    {
        m_theta = 1;
    }


    gsSpaceTimeNorm(const gsField<T> & _field1) 
    : gsNorm<T>(_field1), dfunc2(NULL), f2param(false)
    {
        m_theta = 1;
    }

public:
    
    T compute(bool storeElWise = false)
    {
        this->apply(*this,storeElWise);
        return m_value;
    }

    void setTheta(real_t theta) {m_theta = theta;}

protected:

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        // Setup Quadrature
        const unsigned d = basis.dim();
        gsVector<index_t> numQuadNodes( d );
        for (unsigned i = 0; i < d; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here
        
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_VALUE | NEED_GRAD_TRANSFORM;
    }
    
    // Evaluate on element.
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                  const gsFunction<T>    & func1,
                  const gsFunction<T>    & _func2,
                  gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        func1.deriv_into(quNodes, f1ders);
        // get the gradients to columns

        // Evaluate second function (defined of physical domain)
        geoEval.evaluateAt(quNodes);
        if(dfunc2==NULL)
        {
            // get the gradients to columns
            _func2.deriv_into(geoEval.values(), f2ders);            
        }
        else
            dfunc2->eval_into(geoEval.values(), f2ders);

        // ** Evaluate function v
        //gsMatrix<T> f2val = func2Param ? func2.deriv(quNodes)
        //: func2.eval( geoEval->values() );      
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element, 
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {
        const int d = element.dim();
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Transform the gradients
            // f1pders : N X dim
            // f1pdersSpace, f1pdersTime : N X (dim-1)
            geoEval.transformGradients(k, f1ders, f1pders);
            Rows _f1pdersSpace = f1pders.topRows(d-1);
            Rows _f1pdersTime =  f1pders.bottomRows(1);
            
            Rows f2dersSpace = f2ders.topRows(d-1);
            Rows f2dersTime =  f2ders.bottomRows(1); 
            
            // h - mesh size
            const T h = element.getCellSize();
            const T weight = quWeights[k] *  geoEval.measure(k);
                 
            // f2ders : N X 1
            sum += weight *((_f1pdersSpace - f2dersSpace.col(k)).squaredNorm() + m_theta* h* (_f1pdersTime - f2dersTime.col(k)).squaredNorm() ); //
                           
        }

        accumulated += sum;
        return sum;
    }
    

    inline T takeRoot(const T v) { return math::sqrt(v);}

private:
    real_t m_theta;

    // first derivative of func2:
    const gsFunction<T> * dfunc2; // If this is NULL a numerical approximation will be used

    using gsNorm<T>::m_value;
    using gsNorm<T>::m_elWise;

    gsMatrix<T> f1ders, f2ders;
    gsMatrix<T> f1pders;
    gsMatrix<T> f1pdersTime, f1pdersSpace;
    gsVector<T> unormal;

    bool f2param;// not used yet
};

template <class T>
class gsSpaceTimeSliceNorm : public gsNorm<T>
{
    friend class gsNorm<T>;
public:

    gsSpaceTimeSliceNorm(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             const std::vector<boundaryInterface>& SliceInterfaces,
             const std::vector<patchSide>& InitialSides,
             const std::vector<patchSide>& TopSides,
             bool _f2param = false)
    :gsNorm<T>(_field1,_func2),volumes(_field1,_func2,_f2param), m_sliceInterfaces(SliceInterfaces)
    ,m_initialSides(InitialSides), m_topSides(TopSides),f2param(_f2param)
    {
    }

    gsSpaceTimeSliceNorm(const gsField<T> & _field1,
             const gsFunction<T> & _func2,
             const gsFunction<T> & _dfunc2,
             const std::vector<boundaryInterface>& SliceInterfaces,
             const std::vector<patchSide>& InitialSides,
             const std::vector<patchSide>& TopSides,
             bool _f2param = false)
    : gsNorm<T>(_field1,_func2), volumes(_field1,_func2,_dfunc2,_f2param), m_sliceInterfaces(SliceInterfaces)
    ,m_initialSides(InitialSides), m_topSides(TopSides),f2param(_f2param)
    {
    }


    gsSpaceTimeSliceNorm(const gsField<T> & _field1,const std::vector<boundaryInterface>& SliceInterfaces,
                         const std::vector<patchSide>& InitialSides,
                         const std::vector<patchSide>& TopSides)
    : gsNorm<T>(_field1),volumes(_field1), m_sliceInterfaces(SliceInterfaces)
    ,m_initialSides(InitialSides), m_topSides(TopSides),f2param(false)
    {
    }

    void setTheta(real_t theta) {volumes.setTheta(theta);}
public:

    T compute(bool storeElWise = false)
    {
        m_value = volumes.compute(storeElWise);
        m_value*=m_value; //square the number again.


        gsMultiPatch<T> mp = field1->patches();
        for ( std::vector<boundaryInterface>::const_iterator it=m_sliceInterfaces.begin();it!=m_sliceInterfaces.end();++it) // *it ---> interface
        {
            const gsGeometry<T> & u1  = static_cast<const gsGeometry<T> &>( field1->function(it->first().patch) );
            const gsGeometry<T> & u2  = static_cast<const gsGeometry<T> &>( field1->function(it->second().patch) );

            const bool reverse = u1.basis().numElements(it->first() .side() ) <
                    u2.basis().numElements(it->second().side() ) ;

            const boundaryInterface & iFace =( reverse ? it->getInverse() : *it );
            T curDist = 0;
            if(!reverse)
            {
                curDist= igaDGDistanceJump( mp.patch(iFace.first().patch),
                                            mp.patch(iFace.second().patch),
                                            u1,
                                            u2,
                                            func2->function(0), func2->function(0),
                                            iFace,
                                            f2param);
            }
            else
            {
                curDist= igaDGDistanceJump( mp.patch(iFace.first().patch),
                                            mp.patch(iFace.second().patch),
                                            u2,
                                            u1,
                                            func2->function(0), func2->function(0),
                                            iFace,
                                            f2param);
            }

           m_value += curDist * curDist;
        }

        for ( std::vector<patchSide>::const_iterator it=m_initialSides.begin();it!=m_initialSides.end();++it) // *it ---> interface
        {
            side = it->side();
            this->apply1(*this, storeElWise, it->patch, it->side() );
        }

        for ( std::vector<patchSide>::const_iterator it=m_topSides.begin();it!=m_topSides.end();++it) // *it ---> interface
        {
            side = it->side();
            this->apply1(*this, storeElWise, it->patch, it->side() );
        }

        m_value = takeRoot(m_value);
        return m_value;
    }


protected:

    T igaDGDistanceJump(const gsGeometry<T>& patch1, const gsGeometry<T>& patch2,
                        const gsGeometry<T>& func1,  const gsGeometry<T>& func2, // approximati solution
                        const gsFunction<T>& v1, const gsFunction<T>& v2,	// exact solution
                        const boundaryInterface & bi, // interface
                        bool v_isParam)
    {
        const int d = func1.parDim();
        GISMO_ASSERT ( d == patch1.geoDim(), "Dimension mismatch" );

        typename gsGeometryEvaluator<T>::uPtr geoEval1(getEvaluator(NEED_VALUE | NEED_MEASURE, patch1));

        typename gsGeometryEvaluator<T>::uPtr geoEval2(getEvaluator(NEED_VALUE | NEED_MEASURE, patch2));
        // assuming real-valued function
        typename gsGeometryEvaluator<T>::uPtr funcEval1(getEvaluator(NEED_VALUE, func1));

        typename gsGeometryEvaluator<T>::uPtr funcEval2(getEvaluator(NEED_VALUE, func2));

        const boxSide side1 = bi.first();
        const boxSide side2 = bi.second();

        // "DG method not implemented yet for non-matching interfaces");

        // assumes matching orientation
        // degree of the underlying Gauss rule to use
        gsVector<int> intNodes1 ( func1.basis().dim() );
        const int dir1 = bi.first().direction();
        for (int i = 0; i < dir1; ++i)
            intNodes1[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;
        intNodes1[dir1] = 1;
        for (int i = dir1+1; i < func1.basis().dim(); ++i)
            intNodes1[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;

        gsVector<int> intNodes2 ( func2.basis().dim() );
        const int dir2 = bi.second().direction();
        for (int i = 0; i < dir2; ++i)
            intNodes2[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;
        intNodes2[dir2] = 1;
        for (int i = dir2+1; i < func1.basis().dim(); ++i)
            intNodes2[i] = ( func1.basis().degree(i) + func2.basis().degree(i) + 2 )/ 2 ;


        // Temporaries
        gsVector<T> unormal(d),unormal2(d);

        T sum(0);
        // iterator on grid cells on the "right"
        typename gsDomainIterator<T>::uPtr domIter2= func2.basis().makeDomainIterator(side2);

        const int bSize1      = func1.basis().numElements( bi.first() .side() );
        const int bSize2      = func2.basis().numElements( bi.second().side() );
        const int ratio = bSize1 / bSize2;

        int count = 0;

        gsGaussRule<T> quRule1(intNodes1);
        gsGaussRule<T> quRule2(intNodes2);
        gsMatrix<T> quNodes1, quNodes2;
        gsVector<T> quWeights1, quWeights2;

        // iterate over all boundary grid cells on the "left"
        for (typename gsDomainIterator<T>::uPtr domIter1 = func2.basis().makeDomainIterator(side1);
             domIter1->good(); domIter1->next())
        {
            count++;
            // Compute the quadrature rule on both sides
            quRule1.mapTo(domIter1->lowerCorner(), domIter1->upperCorner(), quNodes1, quWeights1);
            quRule2.mapTo(domIter2->lowerCorner(), domIter2->upperCorner(), quNodes2, quWeights2);

            // compute image of Gauss nodes under geometry mapping as well
            // as Jacobians
            geoEval1->evaluateAt(quNodes1);
            geoEval2->evaluateAt(quNodes2);

            funcEval1->evaluateAt(quNodes1);
            funcEval2->evaluateAt(quNodes2);

            gsMatrix<T> func1_vals = funcEval1->values();// (!) coping
            gsMatrix<T> func2_vals = funcEval2->values();// (!) coping

            // exact solution
            gsMatrix<T> v1_vals = v_isParam ? v1.eval(quNodes1)
                                            : v1.eval( geoEval1->values() );

            gsMatrix<T> v2_vals = v_isParam ? v2.eval(quNodes2)
                                            : v2.eval( geoEval2->values() );

            for (index_t k=0; k!= quWeights1.size(); ++k)
            {
                // Compute the outer normal vector from patch1
                geoEval1->outerNormal(k, side1, unormal);


                // Integral transformation and quadarature weight (patch1)
                // assumed the same on both sides
                const T fff = quWeights1[k] *  unormal.norm();

                const T diff = (func1_vals(0,k) - v1_vals(0,k)) - (func2_vals(0,k) - v2_vals(0,k)) ;
                sum += fff * diff*diff;

            }
            if ( count % ratio == 0 ) // next master element ?
            {
                domIter2->next();
            }
        }
        return takeRoot(sum);
    }

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        // Setup Quadrature
        const unsigned d = basis.dim();
        const int dir = side.direction();
        gsVector<index_t> numQuadNodes( d );
        for (unsigned i = 0; i < d; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE| NEED_VALUE;
    }

    // Evaluate on element.
    inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                         const gsFunction<T>    & _func1,
                         const gsFunction<T>    & _func2,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        _func1.eval_into(quNodes, f1vals);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate second function (defined of physical domain)
        _func2.eval_into(f2param? quNodes : geoEval.values(), f2vals);
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            geoEval.outerNormal(k, side, unormal);
            const T weight = quWeights[k] * unormal.norm();
            sum +=weight * ( f1vals.col(k) - f2vals.col(k) ).squaredNorm();
        }
        accumulated += sum;
        return sum;
    }

    inline T takeRoot(const T v) { return math::sqrt(v);}

private:

    gsSpaceTimeNorm<T> volumes;
    const std::vector<boundaryInterface>& m_sliceInterfaces;
    const std::vector<patchSide>& m_initialSides;
    const std::vector<patchSide>& m_topSides;

    using gsNorm<T>::m_value;
    using gsNorm<T>::m_elWise;

    gsMatrix<T> f1vals, f2vals;
    gsVector<T> unormal;

    using gsNorm<T>::field1;
    using gsNorm<T>::func2;

    bool f2param;
    boxSide side;
};



} // namespace gismo

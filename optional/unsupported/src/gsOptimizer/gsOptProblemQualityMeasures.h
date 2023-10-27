/** @file gsOptProblemQualityMeasures.h

    @brief Provides the code for optimizing quality measures with ipopt

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsOptimizer/gsOptProblem.h>
#include <gsMSplines/gsMappedSpline.h>
#include <gsRecipeAssembler/gsQMOptOperators.h>

namespace gismo
{


template <typename T>
class gsOptProblemQualityMeasures : public gsOptProblem<T>
{
public:

    gsOptProblemQualityMeasures(const gsMappedSpline<2,T> & geo,const gsQualityMeasureWeights & weights) :
        m_geoDim(geo.geoDim()),m_weights(weights),m_compGeom(geo)
    {
        // Shortcut to the input coefficients
        gsMatrix<T> coefs = m_compGeom.coefs();
        coefs.resize(coefs.rows()*m_geoDim,1);

        // fill multiBasis
        gsMappedBasis<2,real_t>& compBasis=m_compGeom.getCompBasis();
        std::vector<gsBasis<real_t> *> basiscontainer;
        for(index_t i = 0;i<compBasis.nPatches();++i)
            basiscontainer.push_back(compBasis.getBase(i).clone());
        m_multiBasis=gsMultiBasis<T>(basiscontainer,compBasis.getTopol());

        gsVector<int> sizes(1);
        sizes(0)=compBasis.size()*2;
        m_mapper=gsDofMapper(sizes);
        std::vector<index_t> boundaryIndices,innerIndices;
        compBasis.boundary(boundaryIndices);
        //compBasis.innerBoundaries(innerIndices);

        boundaryIndices.insert(boundaryIndices.end(), innerIndices.begin(), innerIndices.end());
        sort( boundaryIndices.begin(), boundaryIndices.end() );
        boundaryIndices.erase( unique( boundaryIndices.begin(), boundaryIndices.end() ), boundaryIndices.end() );

        gsVector<unsigned> boundaryIndicesVec(m_geoDim*boundaryIndices.size());
        for(unsigned i = 0;i<boundaryIndices.size();++i)
        {
            boundaryIndicesVec(i)=boundaryIndices[i];
            boundaryIndicesVec(i+boundaryIndices.size())=boundaryIndices[i]+compBasis.size();
        }
        m_mapper.markBoundary(0,boundaryIndicesVec);
        m_mapper.finalize();

        m_numDesignVars  = m_mapper.freeSize();
        m_numConstraints = 0;
        m_numConJacNonZero = 0;

        // Set boundary coefficients equal to the input once and
        // forever
        m_tmpCoefs.resize( m_geoDim*compBasis.size(), 1);
        for (index_t i = 0; i < boundaryIndicesVec.rows(); ++i) // copy to m_tmpCoefs
            m_tmpCoefs( boundaryIndicesVec(i),0 ) = coefs( boundaryIndicesVec(i),0 );

        // Design bounds
        m_desLowerBounds.setConstant(m_numDesignVars, -1.0e19);
        m_desUpperBounds.setConstant(m_numDesignVars,  1.0e19);

        m_conLowerBounds.resize(m_numConstraints);
        m_conUpperBounds.resize(m_numConstraints);

        // Current design -- used by gsOptProblem as starting point
        m_curDesign.resize(m_numDesignVars, 1);
        // view current design as a control point matrix
        for ( int i = 0; i< m_tmpCoefs.rows(); i++ ) // TO DO: more efficient
        {
            // get index of the design control point
            const index_t ii = m_mapper.index(i,0);

            if ( m_mapper.is_free_index(ii) ) // interior node?
                m_curDesign(ii,0) = coefs(i,0);
        }

        //
        m_conJacRows.resize(m_numConJacNonZero);
        m_conJacCols.resize(m_numConJacNonZero);

        m_minVal=std::numeric_limits<real_t>::max();
        m_minCurDesign = m_curDesign;
    }

    ~gsOptProblemQualityMeasures()
    {
    }

public:

    T evalObj( const gsAsConstVector<T> & u ) const
    {
        // return the value of the objective function
        updateTempDesign(u);
        return m_val;
    }

    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
    {
        updateTempDesign(u);
        result=m_deriv;
    }

    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
    {
    }

    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
    {
    }

    void getFullCoefs(gsMatrix<T>& result)
    {
        result=m_tmpCoefs;
        result.resize(m_compGeom.getCompBasis().size(),m_geoDim);
    }

    void getFullCoefsOfMin(gsMatrix<T>& result)
    {
        gsAsConstVector<T> vec(m_minCurDesign.data(),m_minCurDesign.rows());
        updateTempDesign(vec);
        result=m_tmpCoefs;
        result.resize(m_compGeom.getCompBasis().size(),m_geoDim);
    }

protected:

    void updateTempDesign( const gsAsConstVector<T> & u) const // ..
    {
        bool changed = false;
        // Update interior coefficients
        for ( index_t i = 0; i< m_tmpCoefs.rows(); i++ ) // TO DO: more efficient
        {
            const index_t ii = m_mapper.index(i,0);

            if ( m_mapper.is_free_index(ii) ) // interior node?
                if(m_tmpCoefs(i,0)!=u(ii))
                {
                    m_tmpCoefs(i,0) = u(ii);
                    changed = true;
                }
        }

        if(changed)
            calculateValAndDeriv();

        if(m_val<m_minVal)
        {
            std::cout << "updated!\n";
            m_minCurDesign=u;
            m_minVal=m_val;
        }
    }

    void calculateValAndDeriv() const
    {
        gsMatrix<T> newCoefs = m_tmpCoefs;
        newCoefs.resize(m_compGeom.getCompBasis().size(),m_geoDim);
        m_compGeom.setCoefs(newCoefs);
        gsMatrix<T> resultEval,resultDeriv;
        //gsMultiPatch<T> geoMP = m_compGeom.exportToPatches();
        //CalculateMeasureAndDeriv(geoMP,m_multiBasis,m_compGeom.getCompBasis().getMapper(),m_mapper,resultEval,&resultDeriv,m_weights);
        //m_deriv= resultDeriv.col(0);
        m_val  = resultEval(0,0);
        m_deriv= resultDeriv;
    }

private:
    /// Dimension of the ambient space
    const int m_geoDim;
    /// Basis of the input parameterization
    gsMultiBasis<T> m_multiBasis;
    /// Mapper from control points to design variables
    gsDofMapper m_mapper;
    /// Temporary used to translate back designs as control points
    mutable gsMatrix<T> m_tmpCoefs;
    /// weights for the optimization
    const gsQualityMeasureWeights m_weights;

    mutable gsMappedSpline<2,real_t> m_compGeom;

    mutable gsMatrix<T> m_deriv;
    mutable T m_val;

    mutable T m_minVal;
    mutable gsMatrix<T> m_minCurDesign;

private:

    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;
};

} // end namespace gismo

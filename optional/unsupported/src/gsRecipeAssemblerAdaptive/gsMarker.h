/** @file gsMarker.h

    @brief DESCRIPTION

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#pragma once

#include <gsAssembler/gsAdaptiveRefUtils.h>
#include <gsRecipeAssembler/gsPhysicalSpace.h>
#include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>

namespace gismo {


/**
 * @brief The gsCellMarker class
 * Abstract class for per cell marking.
 * Derived classes implement a specific marking strategy by overloading
 * the computeMarkedCells function.
 * The marked cells will be written to m_markedCells.
 */
class gsMarker
{
private:
    bool m_upToDateFlag;
protected:
    gsMatrix<real_t> m_markedCells;
public:
    virtual ~gsMarker() {}

    void reset() { m_upToDateFlag=false; }
    void mark(const gsMatrix<real_t> values) { compute(values); m_upToDateFlag=true; }

    const gsMatrix<real_t>& getMarked() const
    {
        if(!m_upToDateFlag)
            GISMO_ERROR("Marked requested before the call to mark.");
        return m_markedCells;
    }
protected:
    virtual void compute(const gsMatrix<real_t> values)
    {
        m_markedCells.resizeLike(values);
        for(index_t i=0;i<values.rows();++i)
            markRow(values.row(i),m_markedCells.row(i));
    }

    virtual void markRow(gsMatrix<real_t>::constRow errorRow,gsMatrix<real_t>::Row markedRow) = 0;
};

class gsMarkerRelativeThreshold : public gsMarker
{
private:
    real_t m_refParameter;
public:
    gsMarkerRelativeThreshold(real_t refParameter)
        : m_refParameter(refParameter)
    {}
protected:
    void markRow(gsMatrix<real_t>::constRow errorRow,gsMatrix<real_t>::Row markedRow)
    {
        real_t threshold = errorRow.maxCoeff()*m_refParameter;
        // Now just check for each element, whether the local error
        // is above the computed threshold or not, and mark accordingly.
        for( index_t i=0; i < errorRow.cols(); i++)
            markedRow(0,i)= ( errorRow(0,i) >= threshold ) ;
    }
};

class gsMarkerAbsoluteThreshold : public gsMarker
{
private:
    real_t m_threshold;
public:
    gsMarkerAbsoluteThreshold(real_t refParameter)
        : m_threshold(refParameter)
    {}
protected:
    void markRow(gsMatrix<real_t>::constRow errorRow,gsMatrix<real_t>::Row markedRow)
    {
        // Now just check for each element, whether the local error
        // is above the computed threshold or not, and mark accordingly.
        for( index_t i=0; i < errorRow.cols(); i++)
            markedRow(0,i)= ( errorRow(0,i) >= m_threshold ) ;
    }
};

class gsMarkerPercentage : public gsMarker
{
private:
    real_t m_refParameter;
public:
    gsMarkerPercentage(real_t refParameter)
        : m_refParameter(refParameter)
    {}
protected:
    void markRow(gsMatrix<real_t>::constRow errorRow,gsMatrix<real_t>::Row markedRow)
    {
        // Compute the index from which the refinement should start in the sorted vector
      index_t idxRefineStart = cast<real_t,index_t>(math::ceil( m_refParameter * errorRow.cols()) );// if ref is positive refine at least an element
        // ...and just to be sure we are in range:
        if( idxRefineStart >= errorRow.cols() )
        {
            markedRow.setConstant(1);
            return;
        }
        else if( idxRefineStart <= 0 )
        {
            markedRow.setConstant(0);
            return;
        }

        std::vector<real_t> errorVec;
        errorVec.reserve(errorRow.cols());
        for(index_t i =0;i<errorRow.cols();++i)
            errorVec.push_back(errorRow(0,i));
        std::sort(errorVec.begin(), errorVec.end());

        real_t threshold = errorVec[ idxRefineStart ];
        for( index_t i=0; i < errorRow.cols(); i++)
            markedRow(0,i)= ( errorRow(0,i) >= threshold ) ;
    }
};

class gsMarkerFraction : public gsMarker
{
private:
    real_t m_refParameter;
public:
    gsMarkerFraction(real_t refParameter)
        : m_refParameter(refParameter)
    {}
protected:
    void markRow(gsMatrix<real_t>::constRow errorRow,gsMatrix<real_t>::Row markedRow)
    {
        std::vector<real_t> errorVec;
        errorVec.reserve(errorRow.cols());
        for(index_t i =0;i<errorRow.cols();++i)
            errorVec.push_back(errorRow(0,i));
        std::sort(errorVec.begin(), errorVec.end(), std::greater_equal<real_t>());

        real_t totalErrorThreshold=errorRow.sum()*m_refParameter;
        real_t accumulator=0,threshold=0;
        size_t i;
        for( i = 0; i < errorVec.size() && accumulator<totalErrorThreshold; ++i)
            accumulator+=errorVec[i];
        threshold=errorVec[i];

        for( index_t i2=0; i2 < errorRow.cols(); i2++)
            markedRow(0,i2)= ( errorRow(0,i2) >= threshold ) ;
    }
};

}

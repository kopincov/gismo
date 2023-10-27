/** @file gsDivConSolution.h

    @brief Evaluates a function expressed in coefficients of divergence preserving basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsGeometryEvaluator.h>


//----TODO----
//document

namespace gismo
{


template<class T>
class gsDivConSolution : public gsFunction<T>
{

public:

    /// @name Constructors
    /// @{

    /// Default empty constructor
    gsDivConSolution() { }

    /// Constructor which copies the given coefficient matrix \a solutionCoefs.
    gsDivConSolution(const gsMatrix<T> & solutionCoefs, 
                     const gsGeometry<T> & geo, 
                     std::vector< gsBasis<T> *> const & basis) 
    : m_solutionCoeff(solutionCoefs), m_geometry(geo), m_basis(basis)
    {
        componentShifts.resize(targetDim());
        componentShifts[0] = 0;
        for (index_t k = 1; k < targetDim(); ++k)
        {
            componentShifts[k] = componentShifts[k-1] + m_basis[k-1]->size();
        }
    }


    /// @}

public:

    /*/** @name Evaluation functions
      @{
    */

    /// Evaluates the solution into result
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        //JS2: 2D Verified, 3D not Verified(Tested)

        const index_t TarDim = targetDim();
        GISMO_ASSERT(m_geometry.geoDim() == TarDim,
                     "Geometric dimension and target dimention not matching!");

        const index_t numPts = u.cols();
        gsMatrix<T> B;
        gsMatrix<index_t> ind;

        result.setZero( TarDim, numPts );

        typename gsGeometryEvaluator<T>::uPtr geoEval(// to do: member
            getEvaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM, m_geometry));
        geoEval->evaluateAt(u);

        const gsMatrix<T> & Jac = geoEval->jacobians();

        for (index_t comp = 0; comp< TarDim; ++comp)
        {
            m_basis[comp]->eval_into(u, B);
            m_basis[comp]->active_into(u, ind);
            for ( index_t j=0; j< numPts ; j++ ) // for all points (columns of u)
            {
                const T det = geoEval->jacDet(j);
                for ( index_t i = 0; i < ind.rows(); ++i ) // for all non-zero basis functions
                {
                    result.col(j) += Jac.block(0, comp+j*TarDim, TarDim, 1)
                        * m_solutionCoeff(ind(i,j) + componentShifts[comp], 0) * B(i,j)/det;
                }
            }
        }
    }

    /// @}


    /// @name Accessors
    /// @{

    /// \brief Returns the basis.
    std::vector<gsBasis<T> * > basis() {return m_basis;}

    /// Dimension \em n of the physical space
    short_t geoDim() const {return m_geometry.geoDim();}

    /// Dimension \em n of the coefficients (control points)
    short_t coefDim() const { return m_geometry.coefDim(); }

    /// Dimension of the absent physical space (overriding gsFunction::targetDim())
    short_t targetDim() const { return m_basis.size();}

    /// Dimension \em d of the parameter domain.
    virtual short_t parDim() const { return m_basis[0]->dim(); }


    /// Dimension \em d of the parameter domain (overriding gsFunction::domainDim()).
    short_t domainDim() const { return m_basis[0]->dim(); }
    /// @}



protected:

    // The coefficient of the solution field (at one patch), size of basis X targetDim.
    gsMatrix<T> m_solutionCoeff;

    // The geometry, need this for the jacobian.
    const gsGeometry<T>  & m_geometry;//(&)

    // The basis (not owned)
    std::vector< gsBasis<T> * > m_basis;

    // shift index for components in m_solutionCoeff;
    gsVector<index_t> componentShifts;


}; // class gsDivConSolution


}

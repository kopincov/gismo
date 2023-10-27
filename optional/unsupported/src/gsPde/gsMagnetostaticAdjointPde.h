
/** @file gsMagnetostaticAdjointPde.h

    @brief Describes a Poisson PDE.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/


#pragma once

#include <gsCore/gsPiecewiseFunction.h>
#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsPde/gsPdeWithCoeff.h>
#include <gsPde/gsPoissonPde.h>

namespace gismo
{

/** @brief
    A magnetostatic adjoint PDE.

    This class describes the adjoint of a magnetostatic PDE, with an arbitrary right-hand side
    function.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsMagnetostaticAdjointPde : public gsPde<T>, public gsPdeWithCoeff<T>
{
public:
    typedef gsPde<T> Base;

    gsMagnetostaticAdjointPde( ) { }


    /// Constructor
    gsMagnetostaticAdjointPde(const gsMultiPatch<T>         &domain,
                              const gsBoundaryConditions<T> & bc,
                              const gsMultiPatch<T>  & stateSolution,
                              const gsFunctionExpr<T>  & ReferenceFunction,
                              const gsPiecewiseFunction<T>  & alpha)
        : Base(domain,bc), m_stateSol(stateSolution), m_reference(ReferenceFunction), m_alpha(alpha)
    { }

    //const gsPiecewiseFunction<T>* rhs() const {return &m_rhs; }
    //const gsFunction<T>* rhs(np = 0) const {return &m_rhs.piece(np); } ??
    const gsFunction<T> *    getReference()      const { return &m_reference; }

    virtual bool isSymmetric () const { gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return true;}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os<<"Variational adjoint magnetostatic equation   \u222B \u03BD \u2207p \u00B7 \u2207v dx = \u222B f v ds,  with:\n";
        os<<"Reference function Bd = "<< m_reference <<".\n";
        return os;
    }

    virtual gsPde<T>* restrictToPatch(unsigned np) const
    {
        gsBoundaryConditions<T> bc;
        m_boundary_conditions.getConditionsForPatch(np,bc);

        const gsMultiPatch<T> & resStateSol = m_stateSol.piece(np);
        gsPiecewiseFunction<T> newAlpha = m_alpha.piece(np);

        return new gsMagnetostaticAdjointPde<T>(m_domain.patch(np),bc,resStateSol,m_reference,newAlpha);
    }

    virtual void getAlphaValue(unsigned np, const gsMatrix<T>& x, gsMatrix<T>& res) const
    {
        m_alpha.piece(np).eval_into(x,res);
    }

    virtual T getAlphaValue(unsigned np)
    {
        gsMatrix<T> result;
        m_alpha.piece(np).eval_into(gsMatrix<T>::Zero(m_domain[np].parDim(),1),result);
        return result(0,0);
    }

    const gsPiecewiseFunction<T>* getAlpha() const
    {
        return &m_alpha;
    }

    // Not needed???
    const gsFunction<T>* diffusion() const { return &m_alpha.piece(0); }

    virtual gsMultiPatch<T>* getStateSolution(unsigned np = 0)
    {
        return  &m_stateSol;
    }


    virtual T getCoeffForIETI(unsigned np) const {

        //if(np==0)
         //   gsInfo<<"Assume const heterogenous coefficient on each patch\n";

        gsMatrix<T> result;
        m_alpha.piece(np).eval_into(gsMatrix<T>::Zero(m_domain[np].parDim(),1),result);
        return result(0,0);
    }

protected:

    using Base::m_domain;
    using Base::m_boundary_conditions;

    using Base::m_unknownDim;

    gsMultiPatch<T> m_stateSol;
    gsFunctionExpr<T> m_reference;
    gsPiecewiseFunction<T> m_alpha;
    //auto m_rhs;
}; // class gsMagnetostaticAdjointPde

} // namespace gismo


/** @file gsMagnetostaticShapeDerivPde.h

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
class gsMagnetostaticShapeDerivPde : public gsPde<T>, public gsPdeWithCoeff<T>
{
public:
    typedef gsPde<T> Base;

    gsMagnetostaticShapeDerivPde( ) { }


    /// Constructor
    gsMagnetostaticShapeDerivPde(const gsMultiPatch<T>         &domain,
                              const gsBoundaryConditions<T> & bc,
                              const gsPiecewiseFunction<T> & scalingFactor,
                              const gsMultiPatch<T>  & stateSolution,
                              const gsMultiPatch<T> & adjointSolution,
                              const gsPiecewiseFunction<T>  & alpha,
                              const gsPiecewiseFunction<T>  & rhs,
                              const int numRhs = 1)
        : Base(domain,bc), m_scaling(scalingFactor), m_stateSol(stateSolution), m_adjointSol(adjointSolution), m_alpha(alpha), m_rhs(rhs), m_numRhs(numRhs)
    { }

    //const gsPiecewiseFunction<T>* rhs() const {return &m_rhs; }
    //const gsFunction<T>* rhs(np = 0) const {return &m_rhs.piece(np); } ??
    const gsPiecewiseFunction<T> *    rhs()      const { return &m_rhs; }

    virtual bool isSymmetric () const { gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return true;}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os<<"Shape derivative equation   (\u2207J, X) = \u222B S\u2080 \u22C5 X + S\u2081 : \u2202 X dx,  with:\n";
        os<<"Reference function Bd = "<< m_rhs <<".\n";
        return os;
    }

    virtual gsPde<T>* restrictToPatch(unsigned np) const
    {
        gsBoundaryConditions<T> bc;
        m_boundary_conditions.getConditionsForPatch(np,bc);

        gsPiecewiseFunction<T> newScaling = m_scaling.piece(np);
        const gsMultiPatch<T> & resStateSol = m_stateSol.piece(np);
        const gsMultiPatch<T> & resAjointSol = m_adjointSol.piece(np);
        gsPiecewiseFunction<T> newAlpha = m_alpha.piece(np);
        gsPiecewiseFunction<T> newRhs = m_rhs.piece(np);

        return new gsMagnetostaticShapeDerivPde<T>(m_domain.patch(np),bc,newScaling,resStateSol,resAjointSol,newAlpha,newRhs,m_numRhs);
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

    const gsPiecewiseFunction<T>* getScaling() const
    {
        return &m_scaling;
    }

    const gsPiecewiseFunction<T>* getAlpha() const
    {
        return &m_alpha;
    }

    const gsPiecewiseFunction<T>* getRhs() const
    {
        return &m_rhs;
    }

    virtual gsMultiPatch<T>* getStateSolution(unsigned np = 0)
    {
        return  &m_stateSol;
    }

    virtual gsMultiPatch<T>* getAdjointSolution(unsigned np = 0)
    {
        return  &m_adjointSol;
    }

    virtual T getCoeffForIETI(unsigned np) const {

        //if(np==0)
         //   gsInfo<<"Assume const heterogenous coefficient on each patch\n";

        gsMatrix<T> result;
        m_scaling.piece(np).eval_into(gsMatrix<T>::Zero(m_domain[np].parDim(),1),result);
        return result(0,0);
    }

    virtual int numRhs() const
    {
        return m_numRhs;
    }

protected:

    using Base::m_domain;
    using Base::m_boundary_conditions;

    using Base::m_unknownDim;

    gsPiecewiseFunction<T> m_scaling;
    gsMultiPatch<T> m_stateSol;
    gsMultiPatch<T> m_adjointSol;
    gsPiecewiseFunction<T> m_alpha;
    gsPiecewiseFunction<T> m_rhs;
    int m_numRhs;
    //auto m_rhs; // alternative idea to give the rhs expression as an input
}; // class gsMagnetostaticShapeDerivPde

} // namespace gismo

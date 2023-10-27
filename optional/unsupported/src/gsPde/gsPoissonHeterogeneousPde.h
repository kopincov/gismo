
/** @file gsPoissonHeterogeneousPde.h

    @brief Describes a Poisson PDE.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/


#pragma once

#include <gsPde/gsPdeWithCoeff.h>
#include <gsPde/gsPoissonPde.h>
#include <gsCore/gsPiecewiseFunction.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{

/** @brief
    A Poisson PDE.

    This class describes a Poisson PDE, with an arbitrary right-hand side
    function.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsPoissonHeterogeneousPde : public gsPoissonPde<T>, public gsPdeWithCoeff<T>
{
public:
    typedef gsPde<T> Base;

    gsPoissonHeterogeneousPde( ) { }


    /// Constructor
    gsPoissonHeterogeneousPde(const gsMultiPatch<T>         &domain,
                              const gsBoundaryConditions<T> & bc,
                              const gsPiecewiseFunction<T>  & rhs,
                              const gsPiecewiseFunction<T>  & alpha,
                              const gsFunction<T>           *sol = NULL)
        : gsPoissonPde<T>(domain,bc,rhs,sol), m_alpha(alpha)
    { }

    const gsFunction<T> *    rhs()      const { return &m_rhs.piece(0); }

    virtual bool isSymmetric () const { gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return true;}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os<<"Poisson's equation  -\u0394u = f ,  with:\n";
        os<<"Source function f= "<< m_rhs <<".\n";
        return os;
    }

    virtual gsPde<T>* restrictToPatch(unsigned np) const
    {
        gsBoundaryConditions<T> bc;
        m_boundary_conditions.getConditionsForPatch(np,bc);

        gsPiecewiseFunction<T> funs;
        gsConstantFunction<T> con(0,m_domain.parDim());
        for(size_t k=0; k<m_domain.nPatches();k++)
            if(k==np)
                funs.addPiece(m_alpha.piece(k));
            else
                funs.addPiece(con);

        return new gsPoissonHeterogeneousPde<T>(m_domain.patch(np),bc,m_rhs,funs);
    }

    virtual void getAlphaValue(size_t np, const gsMatrix<T>& x, gsMatrix<T>& res) const
    {
        m_alpha.piece(np).eval_into(x,res);
    }

    virtual T getAlphaValue(size_t np)
    {
        gsMatrix<T> result;
        m_alpha.piece(np).eval_into(gsMatrix<T>::Zero(m_domain[np].parDim(),1),result);
        return result(0,0);
    }

    const gsPiecewiseFunction<T>* getAlpha() const
    {
        return &m_alpha;
    }
    
    const gsFunction<T>* diffusion() const { return &m_alpha.piece(0); }
    
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
    using gsPoissonPde<T>::m_rhs;

    gsPiecewiseFunction<T> m_alpha;
}; // class gsPoissonHeterogeneousPde

} // namespace gismo


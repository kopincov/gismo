
/** @file gsMagnetostaticPde.h

    @brief Describes a magnetostatic PDE.

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
    A magnetostatic PDE.

    This class describes a magnetostatic PDE, with an arbitrary right-hand side
    function.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsMagnetostaticPde : public gsPoissonPde<T>, public gsPdeWithCoeff<T>
{
public:
    typedef gsPde<T> Base;

    gsMagnetostaticPde( ) { }


    /// Constructor
    gsMagnetostaticPde(const gsMultiPatch<T>         &domain,
                              const gsBoundaryConditions<T> & bc,
                              const gsPiecewiseFunction<T>  & rhs,
                              const gsPiecewiseFunction<T>  & alpha,
                              const gsFunction<T>           *sol = NULL)
        : gsPoissonPde<T>(domain,bc,rhs,sol), m_alpha(alpha)
    { }

    //const gsPiecewiseFunction<T>* rhs() const {return &m_rhs; }
    //const gsFunction<T>* rhs(np = 0) const {return &m_rhs.piece(np); } ??
    const gsFunction<T> *    rhs()      const { return &m_rhs.piece(0); }

    virtual bool isSymmetric () const { gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return true;}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        //os<<"Magnetostatic equation  -div(\u03BD \uD835\uDFA9 u) = f ,  with:\n";
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
            if(k==0) // for the expression assembler
                funs.addPiece(m_alpha.piece(np));
            else
                funs.addPiece(con);

        gsPiecewiseFunction<T> newRhs;
        newRhs.addPiece(m_rhs.piece(np));
        return new gsMagnetostaticPde<T>(m_domain.patch(np),bc,newRhs,funs);
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

    virtual int numRhs() const
    {
        return 1;
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
}; // class gsMagnetostaticPde

} // namespace gismo


/** @file gsSpaceTimePoissonPde.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    Created on: 2017-06-06
*/

#pragma once

#include <gsPde/gsPoissonHeterogeneousPde.h>

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
class gsSpaceTimePoissonPde : public gsPoissonHeterogeneousPde<T>
{
      typedef gsPde<T> Base;

public:



    gsSpaceTimePoissonPde( ) { }


    /// Constructor
    gsSpaceTimePoissonPde(const gsMultiPatch<T>         &ST_Domain,
                              const gsBoundaryConditions<T> & bc,
                              const gsPiecewiseFunction<T>  & rhs,
                              const gsFunction<T>  & rhs_x,
                              const gsFunction<T>  & rhs_t,
                              const gsPiecewiseFunction<T>  & alpha_x,
                              const gsFunction<T>           *sol = NULL)
        : gsPoissonHeterogeneousPde<T>(ST_Domain,bc,rhs,alpha_x,sol), m_rhs_x(rhs_x), m_rhs_t(rhs_t)
    {

    }

    virtual bool isSymmetric () const { gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return false;}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os<<"Poisson's equation dt u -\u0394u = f ,  with:\n";
        os<<"Source function f= "<< m_rhs <<".\n";
        return os;
    }

    const gsFunction<T> *    rhs_x()      const { return &m_rhs_x; }
    const gsFunction<T> *    rhs_t()      const { return &m_rhs_t; }


    virtual gsPde<T>* restrictToPatch(unsigned np) const
    {
        gsBoundaryConditions<T> bc;
        m_boundary_conditions.getConditionsForPatch(np,bc);

        gsPiecewiseFunction<T> funs;
        gsConstantFunction<T> con(0,m_domain.parDim()-1);
      //  for(index_t k=0; k<m_domain.nPatches();k++)
          //  if(k==np)
                funs.addPiece(m_alpha.piece(np));
          //  else
            //    funs.addPiece(con);

        //TODO: one has to be carfull, when selecting the patch!
        return new gsSpaceTimePoissonPde<T>(m_domain.patch(np),bc,m_rhs,m_rhs_x,m_rhs_t,funs);
    }

    virtual void getAlphaValue(unsigned np, const gsMatrix<T>& x, gsMatrix<T>& res) const
    {
        m_alpha.piece(np).eval_into(x,res);
    }

    virtual T getAlphaValue(unsigned np)
    {
        gsMatrix<T> result;
        m_alpha.piece(np).eval_into(gsMatrix<T>::Zero(m_domain[np].parDim()-1,1),result);
        return result(0,0);
    }


    virtual T getCoeffForIETI(unsigned np) const {

        gsMatrix<T> result;
        m_alpha.piece(np).eval_into(gsMatrix<T>::Zero(m_domain[np].parDim()-1,1),result);
        return result(0,0);
    }



protected:

    using Base::m_domain;
    using Base::m_boundary_conditions;

    using Base::m_unknownDim;
    using gsPoissonHeterogeneousPde<T>::m_rhs;
    using gsPoissonHeterogeneousPde<T>::m_alpha;

    const gsFunction<T>& m_rhs_x;
    const gsFunction<T>& m_rhs_t;


}; // class gsSpaceTimePoissonPDE




} // namespace gismo


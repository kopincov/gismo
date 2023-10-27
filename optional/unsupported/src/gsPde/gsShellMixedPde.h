/** @file gsShellMixedPde.h

    @brief Pde for a mixed formulation of Kirchhoff-Love shell.


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
    Created on: 2017-06-06
*/

#pragma once

#include <gsPde/gsPde.h>
#include <gsCore/gsPiecewiseFunction.h>
#include <gsPde/gsPointLoads.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/** @brief
    Mixed formulation of Kirchhoff-Love shells

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T=real_t>
class gsShellMixedPde : public gsPde<T>
{
    typedef gsPde<T> Base;

public:
    enum corner_type
    {
        sf = 1,
        ff = 2,
        ss = 3,
        none =4,
        other = 5,

    };

    gsShellMixedPde( ) { }

    /// Constructor
    gsShellMixedPde(const gsMultiPatch<T>         &domain,
                    const gsBoundaryConditions<T> & bc,
                    const gsFunction<T>& force,
                    const gsPointLoads<T>& pLoads,
                    const T E,
                    const T nu,
                    const T thickness,
                    const std::vector<std::vector<corner_type> >& corners,
                    const std::vector<std::vector<gsMatrix<> > >& cornerCouplingCoefs,
                    const std::vector<std::vector<int > >& indFreeComp
                    )
    :gsPde<T>(domain,bc), m_force(force), m_pLoads(pLoads), m_E(E), m_nu(nu), m_thickness(thickness), m_corners(corners), m_cornerCouplingCoefs(cornerCouplingCoefs), m_indFreeComp(indFreeComp)
    {
        m_lambda = m_E * m_nu / ( (1.+m_nu)*(1.-2.*m_nu)) ;
        m_mu     = m_E / (2.*(1.+m_nu)) ;
    }

    const gsFunction<T> * force() const { return &m_force;}
    const gsPointLoads<T>& pLoads() const {return m_pLoads;}

    virtual std::ostream &print(std::ostream &os) const
    {
        os<<"ShellMixedPde";
        return os;
    }



protected:

    using Base::m_domain;
    using Base::m_boundary_conditions;
    using Base::m_unknownDim;

    const gsFunction<T>& m_force;
    const gsPointLoads<T>& m_pLoads;

    
public:
    T m_E;
    T m_nu;
    T m_lambda;
    T m_mu;
    T m_thickness;

    std::vector<std::vector<corner_type> > m_corners;
    std::vector<std::vector<gsMatrix<> > > m_cornerCouplingCoefs;
    std::vector<std::vector<int > > m_indFreeComp;


}; // class gsShellMixedPde




} // namespace gismo


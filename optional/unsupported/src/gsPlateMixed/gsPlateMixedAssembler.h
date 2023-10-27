/** @file gsPlateMixedAssembler.h

    @brief Provides assembler for mixed formulation of Kirchhoff-Love plates.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsShellMixedPde.h>

#include <gsCore/gsForwardDeclarations.h>


namespace gismo
{

/** @brief
    Assembles system matrices and right-hand side for a mixed formulation of Kirchhoff-Love plates.
    
    \ingroup Assembler
*/
template <class T=real_t>
class gsPlateMixedAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

    enum phiBcMethod
    {
        nitsche = 1,
        lagrange  = 2,
        nitscheDerivative = 3,
    };

public:

    gsPlateMixedAssembler()
    { }

    /** @brief Main Constructor of the assembler object.

    \param[in] pde: a mixed formulation of KL plates
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] options
    */

    gsPlateMixedAssembler( const gsShellMixedPde<T>      & pde,
                           const gsMultiBasis<T>         & bases,
                           const gsOptionList & options)
    {
        Base::initialize(pde, bases, options);
    }


    gsPlateMixedAssembler( const gsShellMixedPde<T>      & pde,
                           const std::vector< gsMultiBasis<T> >         & bases,
                           const gsOptionList & options)
    {
        Base::initialize(pde, bases, options);
    }

    virtual gsPlateMixedAssembler<T>* clone() const
    {
        return new gsPlateMixedAssembler<T>(*this);
    }
    virtual gsAssembler<T>* create() const
    {
        return new gsPlateMixedAssembler<T>();
    }

    // Refresh routine
    virtual void refresh();

    // Main assembly routine
    virtual void assemble();

    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() might return a lower diagonal
    /// matrix, if we exploit possible symmetry during assembly
    /// (check: m_matrix.symmetry() == true )
    Eigen::SparseSelfAdjointView< typename gsSparseMatrix<T>::Base, Lower> fullMatrix()
    {
        return m_system.matrix().template selfadjointView<Lower>();
    }


    gsMatrix<T> applyVisitorPsiq(index_t patchIndex, boxSide side, const gsVector<T> &upperCorner, gsDomainIterator<T>& element);
    static gsOptionList defaultOptions();



protected:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

public:
    gsMatrix<T> psiq;
    gsMatrix<index_t> activesFreeBoundary;
    gsMatrix<index_t> activesFreeBoundaryMap;
    index_t numActivesFreeBoundary;

};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPlateMixedAssembler.hpp)
#endif

/** @file gsMasorny.h

    @brief Provides assembler for the Masonry problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, Y. Xia
*/

#pragma once

#include <gsAssembler/gsAssemblerBase2.h>
#include <gsAssembler/gsAssemblerOptions.h>

#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals

#include <gsAssembler/gsPoissonAssembler.h>
//#include <gsAssembler/gsAssemblerUtils.h>

#include <gsPde/gsNewtonIterator.h>

#include <gsPde/gsPointLoads.h>


namespace gismo
{

/** @brief
    Implementation of an (multiple right-hand side) Poisson assembler.

    The Poisson equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).
    
    \ingroup Assembler
*/
template <class T>
class gsMasonry_simple : public gsPoissonAssembler<T>
{
public:
    typedef gsPoissonAssembler<T> Base;

public:

    /** @brief
        Constructor of the object.
    */
    gsMasonry_simple( gsMultiPatch<T> const         & patches, // planar domain
               gsBoundaryConditions<T> const & bconditions,
               const gsFunction<T>           & extForce )
    : Base(patches, gsMultiBasis<T>(patches), bconditions, extForce, dirichlet::elimination, iFace::glue)
    {

    }



public:

    void assemble(const gsMultiPatch<T> & curSolution);

    using Base::assemble;

    void  updateSolution(const gsMatrix<T>& solVector, 
                         gsMultiPatch<T>& result) const;

protected:

    // Members from gsAssemblerBase2
    using Base::m_patches;
    using Base::m_bases;
    using Base::m_bConditions;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_matrix;
    using Base::m_rhs;
    using Base::m_dofs;

    // Members from gsPoissonAssembler.h
    using gsPoissonAssembler<T>::m_rhsFun;
};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////

#include <gsSelfSuppSurf/gsMasonry_simple.hpp>
//#ifndef GISMO_BUILD_LIB
//#include GISMO_HPP_HEADER(gsMasonry.hpp)
//#endif

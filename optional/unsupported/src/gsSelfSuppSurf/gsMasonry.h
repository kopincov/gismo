/** @file gsMasorny.h

    @brief Provides assembler for the Masonry problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):  Y. Xia, A. Mantzaflaris.
*/
#pragma once
#include <gsAssembler/gsCDRAssembler.h>
#include <gsPde/gsPointLoads.h>
namespace gismo
{
    /** @brief
    Implementation of an (multiple right-hand side) masonry assembler.

    The masonry equation: 
    \[ - \left[ {{k_1}\frac{{{\partial ^2}z}}{{\partial {x^2}}} 
    + 2{k_2}\frac{{{\partial ^2}z}}{{\partial x\partial y}} 
    + {k_3}\frac{{{\partial ^2}z}}{{\partial {y^2}}}} \right] 
    = {k_4}\sqrt {1 + {z_{,x}}^2 + {z_{,y}}^2} \]

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).
    
    \ingroup Assembler
*/

template <class T>
class gsMasonry: public gsCDRAssembler<T>
{
public:
    typedef gsCDRAssembler<T> Base;

public:

    /** @brief
        Constructor of the object.
    */
    gsMasonry( gsMultiPatch<T> const         & patches, // planar domain
               gsBoundaryConditions<T> const & bconditions,
               const gsFunction<T>           & extForce,
               const gsFunction<T>           & coeff_A, 
               const gsPointLoads<T>         & pLoads = gsPointLoads<T>() )
    : Base(patches, gsMultiBasis<T>(patches), bconditions, extForce, 
           coeff_A, coeffB, coeffC,
           dirichlet::elimination, iFace::glue, stabilizerCDR::none)
    {
        extload_ptr = &extForce;
        coeff_ptr = &coeff_A;
        coeffB = gsConstantFunction<T> (0,0,2);
        coeffC = gsConstantFunction<T> (0,2);

        m_pLoads = pLoads;
    }

public:

    void assemble();

    void assemble(const gsMultiPatch<T> & curSolution);

    void assembleSystem(const gsMultiPatch<T> & curSolution);

    void applyLoads();

    using Base::assemble;

    void  updateSolution(const gsMatrix<T>& solVector, 
                         gsMultiPatch<T>& result) const;

//    void CalGeoHesnMatr(gsMatrix<T> & GeoHesnMatr);


protected:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;

public:
    const gsFunction<T> *extload_ptr;
    const gsFunction<T> *coeff_ptr;
    gsConstantFunction<T> coeffB;
    gsConstantFunction<T> coeffC;

    gsPointLoads<T>  m_pLoads;

    };

}//namespace gismo

#include <gsSelfSuppSurf/gsMasonry.hpp>



//convDiff project

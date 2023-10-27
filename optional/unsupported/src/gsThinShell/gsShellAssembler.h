/** @file gsShellAssembler.h

    @brief Provides system matrices for the elasticity problem on thin shells.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Goyal, A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsAssemblerBase2.h>
#include <gsPde/gsPointLoads.h>

namespace gismo
{

/** @brief Assembles system matrices for thin shell linear and nonlinear elasticity problems.

    \tparam T coefficient type

    \ingroup ThinShell   
*/
template <class T>
class gsShellAssembler : public gsAssemblerBase2<T>
{
public:
    typedef gsAssemblerBase2<T> Base;
    typedef std::vector<std::pair<patchSide,index_t> > clamped_t;
public:

/** @brief Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] thickness is half of the shell thickness
    \param[in] E_modulus is the E-modulus of the chosen material
    \param[in] poissons_ratio is the poisson ratio of the chosen material
    \param[in] bconditions is a gsBoundaryConditions object describing the boundary conditions
    \param[in] surface_force is a gsFunction object describing the 3D surface force
                and/or 2*\em thickness *body_force
    \param[in] clamped contains pairs of patchside objects and x,y or z direction aimed at getting
                zero derivative for that direction
    \param[in] pLoads is a gsPointLoads object describing the points and corresponding loads
    
    \ingroup Assembler
*/
    gsShellAssembler( gsMultiPatch<T> const & patches,
                    // thickness
                    T thickness,
                    // material properties
                    T E_modulus,
                    T poissons_ratio,
                    // Boundary conditions
                    gsBoundaryConditions<T> const & bconditions,
                    // Surface force (in 3D)
                    const gsFunction<T> & surface_force,
                    // z-clamped side
                    const clamped_t & clamped, 
                    const gsPointLoads<T> & pLoads = gsPointLoads<T>() 
                    );
    
/** @brief Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] thickness is half of the shell thickness
    \param[in] rho is the material mass density
    \param[in] E_modulus is the E-modulus of the chosen material
    \param[in] poissons_ratio is the poisson ratio of the chosen material
    \param[in] bconditions is a gsBoundaryConditions object describing the boundary conditions
    \param[in] surface_force is a gsFunction object describing the 3D surface force
                and/or 2*\em thickness *body_force
    \param[in] clamped contains pairs of patchside objects and x,y or z direction aimed at getting
                zero derivative for that direction
    \param[in] pLoads is a gsPointLoads object describing the points and corresponding loads
    
    \ingroup Assembler
*/
    gsShellAssembler( gsMultiPatch<T> const & patches,
                    // thickness
                    T thickness,
                    // material properties
                    T rho,
                    T E_modulus,
                    T poissons_ratio,
                    // Boundary conditions
                    gsBoundaryConditions<T> const & bconditions,
                    // Surface force (in 3D)
                    const gsFunction<T> & surface_force,
                    // z-clamped side
                    const clamped_t & clamped, 
                    const gsPointLoads<T> & pLoads = gsPointLoads<T>()
                    );

public:

    /// Main assembly routine for the linear case or first Newton step in nonlinear case
    /// It assumes that the shell starts from an undeformed initial configuration.
    void assemble();

    /// Main assembly routine for the non-linear case
    void assemble(const gsMultiPatch<T> & deformed);


    /// Reconstruct solution from computed solution vector
    void constructSolution(const gsMatrix<T>& solVector, 
                           gsMultiPatch<T>& result) const;

    /// Reconstruct solution by adding the computed solution vector to the latest solution
    void updateSolution(const gsMatrix<T>& solVector, 
                        gsMultiPatch<T>& result) const;

    void computeDirichletDofs()
    {

        index_t numDirichlet = 0;
        for (int i = 0; i < 3; ++i)
            numDirichlet += m_dofMappers[i].boundarySize();
        m_ddof.setZero(numDirichlet, 1);
        //computeDirichletDofsIntpl();
    }
    
    /// @brief Returns the global mass matrix
    const gsSparseMatrix<T> & massMatrix() const { return m_massMatrix; }
    
protected:

    /// initializes the pde specfic stuff
    void initializePdeSpecific();

    /// Initializes mappings from the patch-dofs to the global system
    /// matrix according to the boundary conditions
    void initMappers(const clamped_t & clamped);

    /// Neumann contributions
    void assembleNeumann();
    
    /// Apply point loads
    void applyLoads();
    
    /// Assemble the mass matrix
    void assembleMassMatrix();

    /// Computes the Dirichlet DoF values by interpolation
    void computeDirichletDofsIntpl();

protected:

    /// @brief half the shell thickness.
    T             m_thickness;
    /// @brief Lame's first parameter.
    T             m_lambda;
    /// @brief Lame's second parameter.
    T             m_mu;
    /// @brief Material density.
    T             m_rho;

    // Boundary conditions
    //gsBoundaryConditions<T> m_bConditions;

    /// Surface force (in 3D)
    const gsFunction<T> * m_surfaceForce;
   
    // Determines how the (fixed) Dirichlet values should be computed
    //dirichlet::values  m_dirValues;

    // Strategy for enforcing Dirichlet DoFs
    //dirichlet::strategy m_dirStrategy;

    const clamped_t & m_clamped;
    
    gsPointLoads<T>  m_pLoads;
    
    /// Mass matrix
    gsSparseMatrix<T> m_massMatrix;

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
};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsShellAssembler.hpp)
#endif

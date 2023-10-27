/** @file gsStokesIterativeSolver.h

    @brief Solver class for the Stokes problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsAssembler/gsPdeAssembler.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsDivConSolution.h>
#include <gsAssembler/gsStokesAssembler.h>
#include <gsSolver/gsMinimalResidual.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsSimplePreconditioners.h>

namespace gismo
{

/** @brief
    TODO Document
*/

/**
  @brief Class to solve stokes problem with iterative methods and precondioners
*/
template<class T = real_t>
class gsStokesIterativeSolver
{
public:
    /// Default empty constructor
    gsStokesIterativeSolver () : m_precon_u(), m_precon_p()  { }

    ///Constructor matrix and rhs
    gsStokesIterativeSolver(gsSparseMatrix<T> matrix, gsMatrix<T> rhs)
    : m_precon_u(), m_precon_p(),m_matrix(matrix), m_rhs(rhs)    {}

    ///Constructor with gsStokesAssembler.
    gsStokesIterativeSolver(gsStokesAssembler<T> * stokesAssembler)
    : m_precon_u(), m_precon_p(), m_stokesAssembler(stokesAssembler)
    {
        initWithStokesAssembler();
    }

public:
    /// @brief Creates the solution fields from the given solution
    /// vector for unknown \a unk on patch \a p.
    gsFunction<T> * reconstructPatchSolution(index_t unk, int p, bool hasMatrixRhs = false) const
    {
        return m_stokesAssembler->reconstructPatchSolution(unk,p,hasMatrixRhs);
    }

    /// @brief Returns the solution of the linear system
    void reconstructSolution(bool hasMatrixRhs = false);

    //void getMGTest(){gsMultiGrid initMultiGridWithStokesAssembler();}


    /*
    /// Set the maximum number of iterations (default: 1000)
    void setMaxIterations(index_t maxIt) {m_solver->setMaxIterations(maxIt);}

    ///Set the tolerance for the error criteria (default: 1e-10)
    void setTolerance(real_t tol) {m_solver->setTolerance(tol);}

    ///The number of iterations needed to reach the error criteria
    int iterations() const { return m_solver->iterations(); }

    ///The error of the iteretive method
    real_t error() const { return m_solver->error(); }
    */

    /// @brief Returns the solution of the linear system
    const gsMatrix<T> &  systemSolution() const { return m_sysSolution; }

    /// solves with Eigen's spareQR direct solver (handles singular matrix)
    void solveWithDirectSolver()
    {
        typename gsSparseSolver<T>::QR  solver;
        solver.analyzePattern( m_matrix );
        solver.factorize     ( m_matrix );
        m_sysSolution = solver.solve( m_rhs );
    }

private:
    // from gsSimplePreconditioners.h
    void setMPreconP(const gsMultiPatch<T> patches, gsMultiBasis<T> basis, index_t numOfSweeps = 1)
    {
        // Assemble the mass matrix for the pressure space
        gsGenericAssembler<T> massConst(patches, basis);
        const gsSparseMatrix<T> & massMatrixBtmp = massConst.assembleMass();

        // Get full matrix (not just lower triangular)
        massMatrixBtmp.cols();
        gsSparseMatrix<T> _mat = massConst.fullMatrix();

        typename gsGaussSeidelOp< gsSparseMatrix<T>, gsGaussSeidel::symmetric >::Ptr sgs =
            gsGaussSeidelOp< gsSparseMatrix<T>, gsGaussSeidel::symmetric >::make( _mat );
        sgs->setNumOfSweeps(numOfSweeps);
        m_precon_p = sgs;

    }

public:
    void solveWithMGandGS(index_t numPreSmo = 2, index_t numPostSmo = 2, index_t numGSSmo = 2)
    {
        switch (m_tarDim)
        {
        case 1:
            m_precon_u = initMultiGrid<1>(numPreSmo, numPostSmo);
            break;
        case 2:
            m_precon_u = initMultiGrid<2>(numPreSmo, numPostSmo);
            break;
        case 3:
            m_precon_u = initMultiGrid<3>(numPreSmo, numPostSmo);
            break;
        default:
            GISMO_ERROR("Dimension is not valid: "<< m_tarDim );
        }

        //m_precon_p = gsGaussSeidelOp<gsSparseMatrix<T>,gsGaussSeidel::symmetric>::make(m_patches, m_basis_p, numGSSmo);
        setMPreconP(m_patches, m_basis_p, numGSSmo);

        gsBlockOp<>::Ptr blockPrec = initBlockPreconditioner();
        gsMinimalResidual<> MinRes(m_matrix, blockPrec);
        gsMatrix<T> x0;
        x0.setZero(m_idofs,1);
        gsInfo <<"Starting iterations..." <<std::endl;
        gsStopwatch time; time.restart();
        MinRes.solve(m_rhs, x0);

        gsInfo <<"done! It took "   << time.stop() << " sec." << std::endl;
        gsInfo <<"residual error: " << MinRes.error() << std::endl;
        gsInfo <<"    iterations: " << MinRes.iterations() << std::endl;

        m_sysSolution = x0;
    }

    void solveWithExactInverse()
    {
        m_precon_u = evalAinv();
        m_precon_p = evalMassInv();

        gsBlockOp<>::Ptr blockPrec = initBlockPreconditioner();
        gsMinimalResidual<> MinRes(m_matrix, blockPrec);
        gsMatrix<T> x0;
        x0.setZero(m_idofs,1);
        gsInfo <<"Starting iterations..." <<std::endl;
        gsStopwatch time; time.restart();
        MinRes.solve(m_rhs, x0);

        gsInfo <<"done! It tok "    << time.stop() << " sec." << std::endl;
        gsInfo <<"residual error: " << MinRes.error() << std::endl;
        gsInfo <<"    iterations: " << MinRes.iterations() << std::endl;

        m_sysSolution = x0;
    }

    /// @brief Find the condition of the canonical preconditioned system
    T conditionNumber(bool removeSingularity = false)
    {
        //Check if inverse is calculated (if not calculate it)
        if (m_inverseA.cols() != m_sz_A && m_inverseA.rows() != m_sz_A)
        {
            evalAinv();
        }
        if (m_inverseMass.cols() != m_sz_B && m_inverseMass.rows() != m_sz_B)
        {
            evalMassInv();
        }
        gsMatrix<T> Binv;
        Binv.setZero(m_idofs,m_idofs);
        Binv.block(0,0,m_sz_A,m_sz_A) = m_inverseA;
        Binv.block(m_sz_A,m_sz_A,m_sz_B,m_sz_B) = m_inverseMass;
        gsMatrix<T> denseSysMatrix = Binv*m_matrix;
        return gsSolverUtils<T>::conditionNumber(denseSysMatrix, removeSingularity);
    }


private:

    void initWithStokesAssembler();

    template<unsigned d>
    gsMultiGridOp<>::Ptr initMultiGrid(index_t numPreSmo = 2, index_t numPostSmo = 2, index_t coarseLvlMG = 50);

    /// Calculates the inverse of the momentum equation and returns a
    /// preconditioner pointer
    typename gsMatrixOp<typename gsMatrix<T>::Base >::Ptr evalAinv()
    {
        //Check if inverse is calculated (if not calculate it)
        if (m_inverseA.cols() != m_sz_A && m_inverseA.rows() != m_sz_A)
        {
            gsMatrix<T> blockA(m_matrix.block(0,0,m_sz_A,m_sz_A));
            m_inverseA = blockA.inverse();
        }
        return makeMatrixOp(m_inverseA);
    }

    /// Calculates the inverse of the mass matrix (for the pressure)
    /// and returns a preconditioner pointer
    typename gsMatrixOp<typename gsMatrix<T>::Base >::Ptr evalMassInv();

    gsBlockOp<>::Ptr initBlockPreconditioner()
    {
        gsBlockOp<>::Ptr blockPrec = gsBlockOp<>::make(2,2);
        blockPrec->addOperator(0,0,m_precon_u);
        blockPrec->addOperator(1,1,m_precon_p);
        return blockPrec;
    }



    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os <<"Stokes problem: -\u03BD\u0394u + \u2207p = f, \n";
        return os;
    }


private:
    //Solver Parameter

    // Kinematic viscosity
    T m_nu;

    //Must use pointers
    gsLinearOperator<>::Ptr m_precon_u;
    gsLinearOperator<>::Ptr m_precon_p;
    //gsMultiGrid precon_u;

    //gsIterativeSolver *m_solver;

    //The stokes assembler object
    gsStokesAssembler<T> * m_stokesAssembler;

    // The multipatch domain (try to remove, flag system)
    gsMultiPatch<T> m_patches;

    // The boundary condition for velocity.
    gsBoundaryConditions<T> m_bcondition;

    // Geometrical transformation type
    ValueTransformationType m_geoTrans;

    dirichlet::strategy m_dirStrategy;
    // Strategy for dealing with patch interface

    // The discretization bases corresponding to \a m_patches and to
    // the number of solution fields that are to be computed
    // m_bases[u][k]: basis for each component of u on patch k
    // (try to remove, flag system)
    std::vector< gsMultiBasis<T> > m_bases_u;

    //Basis for pressure
    // (try to remove, flag system)
    gsMultiBasis<T> m_basis_p;

    // Global matrix
    gsSparseMatrix<T> m_matrix;

    //Storage for the inverse of A and mass matrix (of pressure)
    //if calculated.
    gsMatrix<T> m_inverseA;
    gsMatrix<T> m_inverseMass;

    // Right-hand side ( multiple right hand sides possible )
    gsMatrix<T>     m_rhs;

    // Solution of the linear system
    gsMatrix<T> m_sysSolution;

    // The Dof mapper is used to map patch-local DoFs to the global DoFs
    // One for each unknown.
    std::vector< gsDofMapper *>  m_dofMapper;

    // *** Information ***
    index_t m_idofs;
    index_t m_sz_A;
    index_t m_sz_B;

    index_t m_tarDim;

}; // class gsStokesIterativeSolver

} // namespace gismo


//#ifndef GISMO_BUILD_LIB
//#include GISMO_HPP_HEADER(gsStokesIterativeSolver.hpp)
#include "gsStokesIterativeSolver.hpp"
//#endif

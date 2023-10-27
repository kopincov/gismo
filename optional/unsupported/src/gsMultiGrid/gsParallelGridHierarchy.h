/** @file gsParallelGridHierarchy.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    Created on: 2018-02-23
*/


#pragma once

#include <gsCore/gsConfig.h>


#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>

#include <gsCore/gsDofMapper.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsOptionList.h>
#include <gsIETI/gsParallelOperator.h>

namespace gismo {



template< typename T >
class gsParallelGridHierarchy
{

public:

    /// @brief This function sets up a multigrid hierarchy by uniform refinement
    ///
    /// @param mBasis                    The gsMultiBasis to be refined (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param assemblerOptions          A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
    /// @param levels                    The number of levels
    /// @param numberOfKnotsToBeInserted The number of knots to be inserted, defaulted to 1
    /// @param multiplicityOfKnotsToBeInserted The multiplicity of the knots to be inserted, defaulted to 1
    ///
    /// \ingroup Solver
    static gsParallelGridHierarchy buildByRefinement(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& assemblerOptions,
        gsSortedVector<size_t> myPatches,
        gsMpiComm comm,
        index_t levels,
        index_t numberOfKnotsToBeInserted = 1,
        index_t multiplicityOfKnotsToBeInserted = 1
        );

    /// @brief This function sets up a multigrid hierarchy by uniform refinement
    ///
    /// @param mBasis                    The gsMultiBasis to be refined (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param options                   A gsOptionList defining the necessary infomation
    ///
    /// \ingroup Solver
    static gsParallelGridHierarchy buildByRefinement(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& options,
        gsSortedVector<size_t> myPatches,
        gsMpiComm comm
        )
    {
        return gsParallelGridHierarchy::buildByRefinement(
            give(mBasis),
            boundaryConditions,
            options,
            myPatches,
            comm,
            options.askInt( "Levels", 3 ),
            options.askInt( "NumberOfKnotsToBeInserted", 1 ),
            options.askInt( "MultiplicityOfKnotsToBeInserted", 1 )
        );
    }


    /// @brief This function sets up a grid hierarchy by coarsening
    ///
    /// @param mBasis                    The gsMultiBasis to be coarsened (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param assemblerOptions          A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
    /// @param levels                    The maximum number of levels
    /// @param degreesOfFreedom          Number of dofs in the coarsest grid in the grid hierarchy
    ///
    /// The algorithm terminates if either the number of levels is reached or the number of degrees of freedom
    /// is below the given threshold.
    ///
    static gsParallelGridHierarchy buildByCoarsening(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& assemblerOptions,
        gsSortedVector<size_t> myPatches,
        gsMpiComm comm,
        index_t levels,
        index_t degreesOfFreedom = 0
        );

    /// @brief This function sets up a grid hierarchy by coarsening
    ///
    /// @param mBasis                    The gsMultiBasis to be coarsened (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param options                   A gsOptionList defining the necessary infomation
    ///
    static gsParallelGridHierarchy buildByCoarsening(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& options,
        gsSortedVector<size_t> myPatches,
        gsMpiComm comm
    )
    {
        return gsParallelGridHierarchy::buildByCoarsening(
            give(mBasis),
            boundaryConditions,
            options,
            myPatches,
            comm,
            options.askInt( "Levels", 3 ),
            options.askInt( "DegreesOfFreedom", 0 )
        );
    }

    /// Get the default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt( "DirichletStrategy", "Method for enforcement of Dirichlet BCs [11..14]", 11 );
        opt.addInt( "InterfaceStrategy", "Method of treatment of patch interfaces [0..3]", 1  );
        opt.addInt( "Levels", "Number of levels to be constructed in the grid hierarchy", 3 );
        opt.addInt( "DegreesOfFreedom",   "Number of dofs in the coarsest grid in the grid hierarchy (only buildByCoarsening)", 0 );
        opt.addInt( "NumberOfKnotsToBeInserted", "The number of knots to be inserted (only buildByRefinement)", 1 );
        opt.addInt( "MultiplicityOfKnotsToBeInserted",   "The multiplicity of the knots to be inserted (only buildByRefinement)", 1 );
        opt.addSwitch( "NoSubassembledOperators", "Prevent usage of subassembled operators", false);
        return opt;
    }

    /// Get the stored options
    const gsOptionList& getOptions() const
    { return m_options; }

    /// Reset the object (to save memory)
    void clear()
    {
        //m_boundaryConditions.clear();
        //m_options.clear();
        m_mBases.clear();
        m_transferMatrices.clear();
    }

    /// Get the vector of multi bases (by reference)
    const std::vector< gsMultiBasis<T> >& getMultiBases() const
    { return m_mBases; }
    /// Get the vector of multi bases
    gsParallelGridHierarchy& moveMultiBasesTo( std::vector< gsMultiBasis<T> >& o )
    { o = give(m_mBases); return *this; }

    /// Get the vector of transfer matrices (by reference)
    const std::vector<std::vector< gsSparseMatrix<T, RowMajor> > >& getTransferMatrices() const
    { return m_transferMatrices; }
    /// Get the vector of transfer matrices
    gsParallelGridHierarchy& moveTransferMatricesTo( std::vector<std::vector< gsSparseMatrix<T, RowMajor> > >& o )
    { o = give(m_transferMatrices); return *this; }

    /// Get the vector of transfer operators
    std::pair<std::pair<std::vector<typename gsParallelOperator<T>::Ptr>, std::vector<typename gsParallelOperator<T>::Ptr> >, std::vector<gsParallelGlobalLocalHandler::Ptr> > getRestrictionAndProlongationOperators();

    /// Get the boundary conditions
    const gsBoundaryConditions<T>& getBoundaryConditions() const
    { return m_boundaryConditions; }

    /// Generates the System matrix on the coarser grids
    std::pair<std::vector<typename gsParallelOperator<T>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > generateGalerkinProjection(const std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr >& A, std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr >& A_coarse, std::vector<gsDofMapper>& localMappers_c, gsDofMapper& globalMapper_c) const;
    std::pair<std::vector<typename gsParallelOperator<T>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > generateGalerkinProjection(const std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr >& A) const
    {
        std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr > dummy;
        std::vector<gsDofMapper> localMappers_c;
        gsDofMapper globalMapper_c;
        return generateGalerkinProjection(A,dummy,localMappers_c,globalMapper_c);
    }

    std::pair<std::vector<typename gsParallelOperator<T>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > generateGalerkinProjection(const typename gsMatrixOp<gsSparseMatrix<T> >::Ptr & A, typename gsMatrixOp<gsSparseMatrix<T> >::Ptr  & A_coarse, std::vector<gsDofMapper>& localMappers_c, gsDofMapper& globalMapper_c) const;
    std::pair<std::vector<typename gsParallelOperator<T>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > generateGalerkinProjection(const typename gsMatrixOp<gsSparseMatrix<T> >::Ptr & A) const
    {
        typename gsMatrixOp<gsSparseMatrix<T> >::Ptr dummy;
        std::vector<gsDofMapper> localMappers_c;
        gsDofMapper globalMapper_c;
        return generateGalerkinProjection(A,dummy,localMappers_c,globalMapper_c);
    }

    const std::vector< typename gsPatchSubassembledTopology<T>::Ptr> & getSubassembledTopology() const {return m_subassTopology;}



private:

    gsBoundaryConditions<T> m_boundaryConditions;
    gsOptionList m_options;

    std::vector< gsMultiBasis<T> > m_mBases;
    std::vector< std::vector<gsSparseMatrix<T, RowMajor> > > m_transferMatrices;
    std::vector< typename gsSparseMatrix<T, RowMajor>::Ptr > m_transferMatricesAss;
    std::vector< typename gsPatchSubassembledTopology<T>::Ptr> m_subassTopology;
    std::vector< typename gsConnectionHandler<T>::Ptr > m_connectionHandlers;
    std::vector< gsParallelGlobalLocalHandler::Ptr > m_globLocHandlers;
    gsSortedVector<size_t> m_myPatches;
    gsMpiComm m_comm;

    //Store the mappers on the coarsest grid
    gsDofMapper m_coarseGlobMapper;
    std::vector<gsDofMapper> m_coarseLocMappers;
};


template< typename T>
class gsCoarseSolverAdapter : public gsLinearOperator<T>
{

public:
    /// Shared pointer for gsSolverOp
    typedef memory::shared_ptr<gsCoarseSolverAdapter> Ptr;

    /// Unique pointer for gsSolverOp
    typedef memory::unique_ptr<gsCoarseSolverAdapter> uPtr;

    gsCoarseSolverAdapter(typename gsLinearOperator<T>::Ptr solv,gsParallelGlobalLocalHandler::Ptr handler) : m_solv(solv), m_handler(handler) {}

    static uPtr make(typename gsLinearOperator<T>::Ptr solv,gsParallelGlobalLocalHandler::Ptr handler)
    {
        return memory::make_unique(new gsCoarseSolverAdapter<T>(solv,handler));
    }

    void apply(const gsMatrix<T>& input, gsMatrix<T>& output) const
    {
        m_handler->buildGlobalVector(input,m_tempI);
        m_solv->apply(m_tempI,m_tempO);
        m_handler->extractLocalVector(m_tempO,output);

    }

    index_t rows() const {return m_handler->localSize();}
    index_t cols() const {return m_handler->localSize();}


protected:

    typename gsLinearOperator<T>::Ptr m_solv;
    gsParallelGlobalLocalHandler::Ptr m_handler;

    mutable gsMatrix<T> m_tempI, m_tempO;




};

template<typename T>
typename gsCoarseSolverAdapter<T>::Ptr constructCoarseSolver(
        const std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr >& mat_coarse,
        const gsParallelGlobalLocalHandler::Ptr& handler_coarse,
        const gsSortedVector<size_t>& myPatches,
        const std::vector<gsDofMapper>& patchLocalMappers,
        const gsDofMapper& globalMapper);

template<typename T>
typename gsCoarseSolverAdapter<T>::Ptr constructCoarseSolver(
        const typename gsMatrixOp<gsSparseMatrix<T> >::Ptr & mat_coarse,
        const gsParallelGlobalLocalHandler::Ptr& handler_coarse);





} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include <gsMultiGrid/gsParallelGridHierarchy.hpp>
#endif

#endif

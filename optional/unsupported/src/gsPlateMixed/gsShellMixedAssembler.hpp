/** @file gsShellMixedAssembler.hpp

    @brief Provides assembler implementation for mixed formulation of Kirchhoff-Love plates.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K.Rafetseder
*/

#include <gsPlateMixed/gsShellMixedAssembler.h>

#include <gsPlateMixed/gsVisitorShellMixed.h>
#include <gsPlateMixed/gsVisitorShellMixedFree.h>
#include <gsPlateMixed/gsVisitorShellMixedSimplySupp.h>
#include <gsPlateMixed/gsVisitorShellMixedLineLoad.h>


#include <gsUtils/gsSortedVector.h>

#include <gsCore/gsFuncData.h>



namespace gismo
{

template<class T>
void gsShellMixedAssembler<T>::refresh()
{

    //GISMO_ASSERT(1==m_bases.size(), "We use the same basis for all unknowns");

    // Set mappers and bases for sparse system
    gsVector<index_t> rowInd(12);
    rowInd[0] = 0; // block p uses dofMapper 0
    rowInd[1] = 1; // block phi1 uses dofMapper 1
    rowInd[2] = 2; // block phi2 uses dofMapper 2

    rowInd[3] = 3; // block N11 uses dofMapper 3
    rowInd[4] = 4; // block N22 uses dofMapper 4
    rowInd[5] = 5; // block N12 uses dofMapper 5

    rowInd[6] = 6; // block u1 uses dofMapper 6
    rowInd[7] = 7; // block u2 uses dofMapper 7
    rowInd[8] = 8; // block u3 uses dofMapper 8

    rowInd[9]  = 9;  // block lagrange multiplier lambdaN uses dofMapper 9
    rowInd[10] = 10; // block lagrange multiplier lambdaT uses dofMapper 10
    rowInd[11] = 11; // block lagrange multiplier lambdaMean uses dof Mapper 11


    gsVector<index_t> colvar(12);
    // all blocks use the same basis
    colvar[0] = 0; // p
    colvar[1] = 0; //phi
    colvar[2] = 0;

    colvar[3] = 2; // N // 2
    colvar[4] = 3;      // 3
    colvar[5] = 1;
    //colvar[3] = 0;
    //colvar[4] = 0;
    //colvar[5] = 0;

    colvar[6] = 0; // u
    colvar[7] = 0;
    colvar[8] = 0;

    colvar[9]  = 1; //lambda
    colvar[10] = 1;
    colvar[11] = 0; //alpha


    // 1. Obtain maps from basis functions to matrix columns and rows
    std::vector<gsDofMapper> mappers(12);

    //empty mapper
    gsVector<index_t> patchDofSizes(1);
    patchDofSizes(0)=0;
    gsDofMapper mapperEmpty(patchDofSizes);
    mapperEmpty.setIdentity(1,0);
    mapperEmpty.finalize();

    if(m_options.getSwitch("Mmixed"))
    {
        // p
        m_bases[colvar[0]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[0], 0, 1);
        // phi
        m_bases[colvar[1]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[1], 1, 1);
        m_bases[colvar[2]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[2], 2, 1);
    }
    else
    {
        mappers[0] = mapperEmpty;
        mappers[1] = mapperEmpty;
        mappers[2] = mapperEmpty;
    }


    //N no boundary condition
    if(m_options.getSwitch("Nmixed"))
    {
        m_bases[colvar[3]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[3], 3, 1);
        m_bases[colvar[4]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[4], 4, 1);
        m_bases[colvar[5]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[5], 5, 1);
    }
    else
    {
        mappers[3] = mapperEmpty;
        mappers[4] = mapperEmpty;
        mappers[5] = mapperEmpty;
    }

    // u
    m_bases[colvar[6]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                 (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                 this->pde().bc(), mappers[6], 6, 1);
    m_bases[colvar[7]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                 (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                 this->pde().bc(), mappers[7], 7, 1);
    m_bases[colvar[8]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                 (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                 this->pde().bc(), mappers[8], 8, 1);



    if(m_options.getInt("phiBcMethodF") == lagrange)
    {
        // lagrange multiplier p + (grad phi)_tn
        m_bases[colvar[9]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[9], 9, 0);
        // lagrange multiplier (grad phi)_tt
        m_bases[colvar[10]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[10], 10, 0);


        // auxiliary mapper: only interior Dofs are free
        gsDofMapper mapper_boundary_eliminated;
        gsBoundaryConditions<> BCs_boundary_eliminated;
        for (typename gsBoxTopology::const_biterator it = m_pde_ptr->patches().bBegin(); it!= m_pde_ptr->patches().bEnd(); ++it)
        {
            patchSide ps = *it;
            BCs_boundary_eliminated.addCondition(ps.patch, ps.side(), condition_type::dirichlet, 0, 0 );
        }

        m_bases[colvar[9]].getMapper(dirichlet::elimination,
                                     iFace::glue,
                                     BCs_boundary_eliminated, mapper_boundary_eliminated, 0, 1);

        // eliminate interior Dofs
        for(index_t i=0; i<mapper_boundary_eliminated.freeSize(); i++)
        {
            std::vector<std::pair<index_t,index_t> > result;
            mapper_boundary_eliminated.preImage(i,result); // Todo: improve perfomance
            mappers[9].eliminateDof(result[0].second, result[0].first);
            mappers[10].eliminateDof(result[0].second, result[0].first);
        }

        mappers[9].finalize();
        mappers[10].finalize();
    }
    else
    {
        mappers[9] = mapperEmpty;
        mappers[10] = mapperEmpty;
    }

    // langrange multiplier mean lambdaT = 0
    if(m_options.getSwitch("lambdaTMean0"))
    {
        int nLambdaTMean0 = 0; // number of Lagrangemult for mean

        gsVector<index_t> patchDofSizes_temp(m_pde_ptr->patches().nPatches());
        patchDofSizes_temp.setOnes(m_pde_ptr->patches().nPatches());
        patchDofSizes_temp *= nLambdaTMean0;

        gsDofMapper mapper(patchDofSizes_temp);
        mapper.setIdentity(m_pde_ptr->patches().nPatches(), nLambdaTMean0);
        mapper.finalize();
        mappers[11] = mapper;
    }
    else
    {
        mappers[11] = mapperEmpty;
    }



    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(mappers, rowInd, rowInd, colvar); // rows and cols use the same mappers

}


template<class T>
void gsShellMixedAssembler<T>::assemble()
{
    GISMO_ASSERT(m_system.initialized(),
                 "Sparse system is not initialized, call initialize() or refresh()");

    // Reserve memory
    const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0], 2, 1, .33); // uses basis[0] of patch 0 (multipatch) !!
    m_system.reserve(9*nz, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    if(m_options.getSwitch("Mmixed"))
    {
        Base::computeDirichletDofs(0);
        Base::computeDirichletDofs(1);
        Base::computeDirichletDofs(2);
    }

    Base::computeDirichletDofs(6);
    Base::computeDirichletDofs(7);
    Base::computeDirichletDofs(8);

    if(m_options.getInt("phiBcMethodF") == lagrange)
    {
        Base::computeDirichletDofs(9);
        Base::computeDirichletDofs(10);
    }


    // Clean the sparse system
    //m_system.setZero(); // !! does not keep allocated memory !!


    // Assemble volume integrals
    std::cout<<"volume integrals"<<std::endl;
    gsVisitorShellMixed<T> visitor(*m_pde_ptr, *this);
    Base::template push<gsVisitorShellMixed<T> >(visitor);
    //Base::template push<gsVisitorShellMixed<T> >();


    // Assemble boundary integrals
    std::cout<<"boundary integrals"<<std::endl;
    if(m_options.getSwitch("Mmixed"))
    {
        pushWithAssembler<gsVisitorShellMixedSimplySupp<T> >(m_pde_ptr->bc().container("SimplySupp") );
        pushWithAssembler<gsVisitorShellMixedFree<T> >(m_pde_ptr->bc().container("Free") );
        pushWithAssembler<gsVisitorShellMixedLineLoad<T> >(m_pde_ptr->bc().container("LineLoad") );

        //Base::template push<gsVisitorShellMixedSimplySupp<T> >(m_pde_ptr->bc().container("SimplySupp") );
        //Base::template push<gsVisitorShellMixedFree<T> >(m_pde_ptr->bc().container("Free") );
        //Base::template push<gsVisitorShellMixedLineLoad<T> >(m_pde_ptr->bc().container("LineLoad") );
    }

    // Apply point loads
    applyPointLoads();

    // Assembly is done, compress the matrix
    Base::finalize();
}



template<class T>
gsOptionList gsShellMixedAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();

    opt.addInt("phiBcMethodSs", "phi problem method for BC on ss part: nitsche 1, Lagrange 2", 1);
    opt.addInt("phiBcMethodF", "phi problem method for BC on f part: nitsche 1, Lagrange 2", 1);
    opt.addSwitch("wHomogenPsiq", "w problem homogenization with psi[q]", 0);
    opt.addSwitch("wPenaltyPsiq", "w problem penalty with p(phi - psi[p], psi[q])", 0);
    opt.addSwitch("solveWholeSystem", "whole system is solved at once", 0);
    opt.addSwitch("pBoundary0", "", 0);
    opt.addSwitch("Mmixed", "auxiliary variable M", 0);
    opt.addSwitch("Nmixed", "auxiliary variable N", 0);
    opt.addSwitch("lambdaTMean0", "", 0);

    return opt;

}

template<class T>
void gsShellMixedAssembler<T>::applyPointLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<index_t> actives, globalActives;
    gsMatrix<T,3,3>    F;
    gsVector<T> normal;
    gsMatrix<T,3,1>    value;
    const gsShellMixedPde<T>& shellPde_ptr= dynamic_cast<const gsShellMixedPde<T>& >(*m_pde_ptr);

    const gsPointLoads<T>& pLoads = shellPde_ptr.pLoads();

    unsigned evFlags = NEED_JACOBIAN | NEED_NORMAL;

    for (size_t i = 0; i< pLoads.numLoads(); ++i )
    {
        // compute contravariant components of forceVals in the covariant basis
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_pde_ptr->patches()[pLoads[i].patch]));
        geoEval->evaluateAt(pLoads[i].point);
        geoEval->normal(0,normal);
        normal.normalize();
        F.leftCols(2) = geoEval->jacobian(0);
        F.col(2)      = normal;
        value = F.inverse() * pLoads[i].value;

        if ( pLoads[i].parametric )
        {
            m_bases.front().basis(pLoads[i].patch).active_into( pLoads[i].point, actives );
            m_bases.front().basis(pLoads[i].patch).eval_into  ( pLoads[i].point, bVals);
        }
        else
        {
            gsWarn<< "Point loads parametric for now.\n";
        }

        // translate patch-local indices to global dof indices
        for (size_t j = 0; j< 3; ++j)
        {
            if (value[j] != 0.0)
            {
                globalActives.setZero(actives.rows(), 1);

                for (index_t k=0; k < actives.rows(); ++k)
                {
                    m_system.mapToGlobalRowIndex(actives(k,0), pLoads[i].patch, globalActives(k,0), j+6); // u1, u2, u3 row index 6,7,8
                    m_system.rhs()(globalActives(k,0), 0) += -bVals(k,0) * value[j];
                }
            }
        }
    }
}

template <class T>
template<class BElementVisitor>
void gsShellMixedAssembler<T>::pushWithAssembler(const bcContainer & BCs)
{
    for (typename bcContainer::const_iterator it
         = BCs.begin(); it!= BCs.end(); ++it)
    {
        BElementVisitor visitor(*m_pde_ptr, *it, *this);
        //Assemble (fill m_matrix and m_rhs) contribution from this BC
        Base::apply(visitor, it->patch(), it->side());
    }
}



}// namespace gismo

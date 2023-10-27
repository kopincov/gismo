/** @file gsPlateMixedAssembler.hpp

    @brief Provides assembler implementation for mixed formulation of Kirchhoff-Love plates.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K.Rafetseder
*/


#include <gsPlateMixed/gsVisitorPlateMixed.h>
#include <gsPlateMixed/gsVisitorPlateMixedFree.h>
#include <gsPlateMixed/gsVisitorPlateMixedSimplySupp.h>
#include <gsPlateMixed/gsVisitorPlateMixedClamped.h>
#include <gsPlateMixed/gsVisitorPlateMixedPsiq.h>

#include <gsUtils/gsSortedVector.h>



namespace gismo
{

template<class T>
void gsPlateMixedAssembler<T>::refresh()
{

    //GISMO_ASSERT(1==m_bases.size(), "We use the same basis for all unknowns");

    // Set mappers and bases for sparse system
    gsVector<index_t> rowInd(7);
    rowInd[0] = 0; // block p uses dofMapper 0
    rowInd[1] = 1; // block phi1 uses dofMapper 1
    rowInd[2] = 2; // block phi2 uses dofMapper 2
    rowInd[3] = 3; // block w uses dofMapper 3
    rowInd[4] = 4; // block lagrange multiplier lambdaN uses dofMapper 4
    rowInd[5] = 5; // block lagrange multiplier lambdaT uses dofMapper 5
    rowInd[6] = 6; // block lagrange multiplier lambdaMean uses dof Mapper 6



    gsVector<index_t> colvar(7);
    // all blocks use the same basis
    colvar[0] = 0;
    colvar[1] = 0;
    colvar[2] = 0;
    colvar[3] = 0;
    colvar[4] = 1; // 1
    colvar[5] = 1; // 1
    colvar[6] = 0;


    // 1. Obtain maps from basis functions to matrix columns and rows
    std::vector<gsDofMapper> mappers(7);

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
    // w
    m_bases[colvar[3]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                 (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                 this->pde().bc(), mappers[3], 3, 1);

    if(m_options.getInt("phiBcMethodF") == lagrange)
    {
        // lagrange multiplier p + (grad phi)_tn
        m_bases[colvar[4]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[4], 4, 0);
        // lagrange multiplier (grad phi)_tt
        m_bases[colvar[5]].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                                     (iFace    ::strategy)m_options.getInt("InterfaceStrategy"),
                                     this->pde().bc(), mappers[5], 5, 0);


        // auxiliary mapper: only interior Dofs are free
        gsDofMapper mapper_boundary_eliminated;
        gsBoundaryConditions<> BCs_boundary_eliminated;
        for (typename gsBoxTopology::const_biterator it = m_pde_ptr->patches().bBegin(); it!= m_pde_ptr->patches().bEnd(); ++it)
        {
            patchSide ps = *it;
            BCs_boundary_eliminated.addCondition(ps.patch, ps.side(), condition_type::dirichlet, 0, 0 );
        }

        m_bases[colvar[4]].getMapper(dirichlet::elimination,
                                     iFace::glue,
                                     BCs_boundary_eliminated, mapper_boundary_eliminated, 0, 1);

        // eliminate interior Dofs
        for(index_t i=0; i<mapper_boundary_eliminated.freeSize(); i++)
        {
            std::vector<std::pair<index_t,index_t> > result;
            mapper_boundary_eliminated.preImage(i,result); // Todo: improve perfomance
            mappers[4].eliminateDof(result[0].second, result[0].first);
            mappers[5].eliminateDof(result[0].second, result[0].first);
        }

        //mappers[4].finalize();
        //mappers[5].finalize();

    }
    else
    {
        gsVector<index_t> patchDofSizes(1);
        patchDofSizes(0)=0;
        gsDofMapper mapperEmpty(patchDofSizes);
        mapperEmpty.setIdentity(1,0);

        mappers[4] = mapperEmpty;
        mappers[5] = mapperEmpty;
        mappers[6] = mapperEmpty;
    }

    // langrange multiplier mean lambdaT = 0
    if(m_options.getSwitch("lambdaTMean0"))
    {
        gsVector<index_t> patchDofSizes(1);
        patchDofSizes(0)=1;
        gsDofMapper mapper_6(patchDofSizes);
        mapper_6.setIdentity(1,1);
        mappers[6] = mapper_6;
    }
    else
    {
        gsVector<index_t> patchDofSizes(1);
        patchDofSizes(0)=0;
        gsDofMapper mapper_6(patchDofSizes);
        mapper_6.setIdentity(1,0);
        mappers[6] = mapper_6;
    }

    mappers[4].finalize();
    mappers[5].finalize();
    mappers[6].finalize();


    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(mappers,rowInd,rowInd,colvar); // rows and cols use the same mappers


    // Find actives on free boundary
    psiq.setZero(m_system.colMapper(0).freeSize(),2);
    gsSortedVector<index_t> activesFreeCorners;

    typename gsBoundaryConditions<T>::bcContainer freeSides = m_pde_ptr->bc().robinSides();
    numActivesFreeBoundary = 0;
    for (typename gsBoundaryConditions<T>::bcContainer::const_iterator it= freeSides.begin(); it!= freeSides.end(); ++it)
    {
        gsMatrix<index_t> bPtr = m_bases[m_system.colBasis(0)][it->patch()].boundary(it->side());
        numActivesFreeBoundary += (bPtr).rows();

        if(it!=freeSides.begin())
            numActivesFreeBoundary -=1;

    }

    activesFreeBoundary.resize(numActivesFreeBoundary, 1);
    gsMatrix<index_t> activesFreeBoundary_temp;
    int curRowInd=0;
    for (typename gsBoundaryConditions<T>::bcContainer::const_iterator it= freeSides.begin(); it!= freeSides.end(); ++it)
    {
        gsMatrix<index_t> bPtr = m_bases[m_system.colBasis(0)][it->patch()].boundary(it->side());
        activesFreeBoundary_temp = bPtr;

        activesFreeBoundary.block(curRowInd, 0, activesFreeBoundary_temp.rows()-2, 1) = activesFreeBoundary_temp.block(1, 0, activesFreeBoundary_temp.rows()-2, 1);
        curRowInd += activesFreeBoundary_temp.rows()-2;
        // corners
        activesFreeCorners.push_sorted_unique(activesFreeBoundary_temp(0,0));
        activesFreeCorners.push_sorted_unique(activesFreeBoundary_temp(activesFreeBoundary_temp.rows()-1,0));
    }

    for(size_t i= 0; i!= activesFreeCorners.uniqueSize(); i++)
    {
        activesFreeBoundary(curRowInd, 0)= activesFreeCorners[i];
        curRowInd++;
    }

    m_system.mapColIndices(activesFreeBoundary, 0, activesFreeBoundaryMap, 0); // TODO: fix multipatch




    // reserve memory

    gsMatrix<index_t> bPtr = m_bases[m_system.colBasis(1)][0].allBoundary(); // use phi basis, because higher degree  // TODO: fix multipatch
    index_t  numActivesBoundary = (bPtr).rows();

    const index_t nz = gsAssemblerOptions::numColNz(m_bases[m_system.colBasis(1)][0], 2, 1, .33); // use phi basis, because higher degree  // TODO: fix multipatch

    // old version (same for all columns)
    //m_system.reserve(4*nz + numActivesFreeBoundary + 2*(m_bases[0][0].maxDegree()+1)*numActivesBoundary, this->pde().numRhs()); // multipatch: basis of patch 0 is considered

    // seperately per column
    gsVector<int> nonZerosPerCol(m_system.cols());

    // bilinear forms
    for(index_t i=0; i<m_system.cols();i++)
    {
        if(i< m_system.colMapper(0).freeSize())
            nonZerosPerCol(i) = 4*nz;
        else if(i< m_system.colMapper(0).freeSize() + 2*m_system.colMapper(1).freeSize())
            nonZerosPerCol(i) = 3*nz;
        else
            nonZerosPerCol(i) = 2*nz; // 1
    }

    nonZerosPerCol(m_system.cols()-1) = m_bases[m_system.colBasis(1)][0].size();


    // psiq terme

    if(m_options.getInt("phiBcMethodF") == nitsche)
    {

        gsMatrix<index_t> activesBoundary;
        for (typename gsBoxTopology::const_biterator it = m_pde_ptr->patches().bBegin(); it!= m_pde_ptr->patches().bEnd(); ++it)
        {
            // loop over column blocks 0, 1, 2
            for(index_t c=0; c<=2; c++)
            {
                // offSet 0,...,p+1 needed to find all functions with support on the elements corresponding to the boundary
                for(int offSet = 0; offSet<=m_bases[m_system.colBasis(c)][it->patch].maxDegree()+1; offSet++) // use phi basis, because fine
                {
                    activesBoundary = m_bases[m_system.colBasis(c)][it->patch].boundaryOffset(it->side(), offSet);

                    for(index_t i=0; i<activesBoundary.rows(); i++)
                    {
                        index_t indMap;
                        m_system.mapToGlobalColIndex(activesBoundary(i,0), it->patch, indMap, c);
                        nonZerosPerCol(indMap) += numActivesFreeBoundary;
                    }
                }
            }
        }


        // psip terme
        for(index_t i=0; i<activesFreeBoundaryMap.rows(); i++)
        {
            nonZerosPerCol(activesFreeBoundaryMap(i,0)) += 3*(m_bases[m_system.colBasis(1)][0].maxDegree()+1)*numActivesBoundary; // use phi basis, because higher degree // TODO: fix multipatch
        }
    }

    m_system.matrix().reserve(nonZerosPerCol);
    m_system.rhs().setZero(m_system.cols(), this->pde().numRhs());



    // m_ddof.resize(tarDim+1); !!
}


template<class T>
void gsPlateMixedAssembler<T>::assemble()
{
    GISMO_ASSERT(m_system.initialized(),
                 "Sparse system is not initialized, call initialize() or refresh()");

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    Base::computeDirichletDofs(0);
    Base::computeDirichletDofs(1);
    Base::computeDirichletDofs(2);
    Base::computeDirichletDofs(3);
    if(m_options.getInt("phiBcMethodF") == lagrange)
    {
        Base::computeDirichletDofs(4);
        Base::computeDirichletDofs(5);
    }
    if(m_options.getSwitch("lambdaTMean0"))
    {

    }

    // Clean the sparse system
    //m_system.setZero(); // !! does not keep allocated memory !!


    // Assemble volume integrals
    std::cout<<"volume integrals"<<std::endl;
    //Base::template push<gsVisitorPlateMixed<T> >();
    for (size_t np=0; np < m_pde_ptr->domain().nPatches(); ++np )
    {
        gsVisitorPlateMixed<T> visitor(*m_pde_ptr, *this);
        Base::apply(visitor, np);
    }



    // Assemble boundary integrals
    std::cout<<"boundary integrals"<<std::endl;
    typename gsBoundaryConditions<T>::bcContainer BCs;
    for (typename gsBoxTopology::const_biterator it = m_pde_ptr->patches().bBegin(); it!= m_pde_ptr->patches().bEnd(); ++it)
    {
        std::cout<<"boundary"<<std::endl;
        std::cout<<*it<<std::endl;
        m_pde_ptr->bc().getConditionFromSide(*it, BCs);

        for (typename gsBoundaryConditions<T>::const_iterator BC = BCs.begin(); BC!= BCs.end(); ++BC)
        {
            if(BC->type() == condition_type::neumann)
            {
                gsVisitorPlateMixedSimplySupp<T> visitor(*m_pde_ptr, *BC, *this);
                Base::apply(visitor, BC->patch(), BC->side());
            }

            else if(BC->type() == condition_type::robin)
            {
                gsVisitorPlateMixedFree<T> visitor(*m_pde_ptr, *BC, *this);
                Base::apply(visitor,BC->patch(), BC->side());
            }

        }
    }


    // Clamped boundary
    BCs = m_pde_ptr->bc().container("Clamped");
    for (typename gsBoundaryConditions<T>::const_iterator BC= BCs.begin(); BC!= BCs.end(); ++BC)
    {
        gsVisitorPlateMixedClamped<T> visitor(*m_pde_ptr, *BC, *this);
        Base::apply(visitor,BC->patch(), BC->side());
    }



    /*
     // If requested, enforce Dirichlet boundary conditions by Nitsche's method
     if ( m_options.dirStrategy == dirichlet::nitsche )
         Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());
     // If requested, enforce Dirichlet boundary conditions by diagonal penalization
     else if ( m_options.dirStrategy == dirichlet::penalize )
         Base::penalizeDirichletDofs();
    */

    if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
        gsWarn <<"DG option is ignored.\n";
    
    // Assembly is done, compress the matrix
    Base::finalize();
}


template <class T>
gsMatrix<T> gsPlateMixedAssembler<T>::applyVisitorPsiq(index_t patchIndex,
                                                       boxSide side, const gsVector<T>& upperCorner, gsDomainIterator<T>& element)
{

    const gsBasisRefs<T> bases(m_bases, patchIndex);

    gsQuadRule<T> QuRule ; // Quadrature rule
    gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
    gsVector<T> quWeights; // Temp variable for mapped weights
    unsigned evFlags(0);

    gsVisitorPlateMixedPsiq<T> visitor(side);

    // Initialize reference quadrature rule and visitor data
    visitor.initialize(bases, patchIndex, m_options, QuRule, evFlags);

    //fixme: gsMapData<T> mapData;
    // Initialize geometry evaluator
    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_pde_ptr->patches()[patchIndex]));

    // Flip corners if necessary (s.t. counter-clockwise)
    gsVector<T> element_lowerCorner(element.lowerCorner());
    gsVector<T> element_lowerCorner_tmp(element.lowerCorner());
    gsVector<T> element_upperCorner(element.upperCorner());
    if((sideOrientation(side) * geoEval->orientation()) == -1)
    {
        if(side.direction() == 0)
        {
            element_lowerCorner[1] = 1 - element_upperCorner[1];
            element_upperCorner[1] = 1 - element_lowerCorner_tmp[1];
        }
        else
        {
            element_lowerCorner[0] = 1 - element_upperCorner[0];
            element_upperCorner[0] = 1 - element_lowerCorner_tmp[0];
        }
    }


    QuRule.mapTo( element_lowerCorner, upperCorner, quNodes, quWeights );

    // Perform required evaluations on the quadrature nodes
    visitor.evaluate(bases, /* *domIt,*/ *geoEval, quNodes);

    // Assemble on element
    visitor.assemble(element, *geoEval, quWeights);

    // Push to global matrix and right-hand side vector
    if(element_upperCorner == upperCorner)
    {
        visitor.localToGlobal(patchIndex, m_system, psiq);
    }

    return visitor.getLocalMat();
}


template<class T>
gsOptionList gsPlateMixedAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();

    opt.addInt("phiBcMethodSs", "phi problem method for BC on ss part: nitsche 1, lagrange 2, nitsche derivative 3", 1);
    opt.addInt("phiBcMethodF", "phi problem method for BC on f part: nitsche 1, lagrange 2, nitsche derivative 3", 1);
    opt.addSwitch("wHomogenPsiq", "w problem homogenization with psi[q]", 0);
    opt.addSwitch("wPenaltyPsiq", "w problem penalty with p(phi - psi[p], psi[q])", 0);
    opt.addSwitch("solveWholeSystem", "whole system is solved at once", 0);
    opt.addSwitch("pBoundary0", "", 0);
    opt.addSwitch("lambdaTMean0", "", 0);

    return opt;

}



}// namespace gismo

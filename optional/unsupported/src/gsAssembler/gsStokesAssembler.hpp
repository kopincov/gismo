/** @file gsStokesAssembler.hpp

    @brief Assembler and solver for the Stokes problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#pragma once

#include <gsCore/gsDebug.h>
#include <gsCore/gsBasisEvaluator.h>

#include <gsUtils/gsBVProblem.h>
#include <gsAssembler/gsAssemblerUtils.h>
#include <gsRecipeAssembler/gsLocalToGlobal.h>
#include <gsMultiGrid/gsGridHierarchy.h>

#include <gsCore/gsGeometryEvaluator.h>


namespace gismo {



template<class T>
void gsStokesAssembler<T>::initDofMapper(iFace::strategy interfaceStrategy,
                                      dirichlet::strategy dirStrategy, bool hasMatrixRhs)
{
    // Setup the degrees of freedom, only implemented for single patch

    // Clean up
    freeAll(m_dofMapper);

    //One for velocity and one for pressure
    m_dofMapper.resize(2);
    const bool match     = (interfaceStrategy == iFace::glue);
    index_t TarDim = m_unknownDim[0];

    //Make a vector fo multibasis with only velocity components.
    std::vector<const gsMultiBasis<T> *> bases_u_ptr(TarDim);
    for (index_t k = 0; k < TarDim; ++k)
    {
        bases_u_ptr[k] = &m_bases[k];

    }
    if (m_geoTrans == DIV_CONFORMING && dirStrategy == dirichlet::eliminatNormal)
    {
        m_dofMapper[0] =  makeDofMapperForVelocityDivConforming(bases_u_ptr, *m_bconditions[0]);
    }
    else
        m_dofMapper[0] =  makeVectorValuedDofMapper(bases_u_ptr, *m_bconditions[0], dirStrategy);
    m_dofMapper[1] = new gsDofMapper;
    m_bases[TarDim].getMapper(match, *m_dofMapper[1]);

    // Set internal and total number of DoFs
    m_idofs = m_dofMapper[0]->freeSize() + m_dofMapper[1]->freeSize();
    m_dofs  = m_dofMapper[0]->size() + m_dofMapper[1]->size();

}

template<class T>
void gsStokesAssembler<T>::initialize()
{
    //Tells the user that elimination strategy does not work in combination with dg.
    GISMO_ASSERT(!(m_interfaceStrategy == iFace::dg && 
                   m_dirStrategy == dirichlet::elimination),
                 "Must use Nitsche method for Dirichlet BC when using discontinuous Galerkin for patch interface");

    m_fixedDofs.resize(2); //One unknown u (velocity) and unknown p (pressure)
    m_fixedDofs[0].resize(0, m_unknownDim[0]);
    index_t tarDim = m_unknownDim[0];

    //Check if each velocity component has the same basis (i.e => A1 = A2 = A3)
    if ((int) m_bases.size() == m_unknownDim.size())
    {
        //Save and remove the bases for the pressure (added at the end)
        gsMultiBasis<T> pressure_bases = m_bases[m_bases.size()-1];
        m_bases.pop_back();
        //Add basis for velocity
        for (index_t k = 0; k < tarDim-1; ++k)
        {
            m_bases.push_back(gsMultiBasis<T>(m_bases[0]));
        }
        m_bases.push_back(pressure_bases);
    }
    else if ((int) m_bases.size() == m_unknownDim.sum())
    {}
    else
        GISMO_ERROR("Number of basis does not match to number of unknowns!");

    //initDofMapper needs to be called before initAssembler
//    GISMO_ASSERT (m_patches.nPatches() == 1, "currently one single patch is working becouse of upgrading of the DOF mapper");
    if (m_geoTrans == DIV_CONFORMING)
    {
        initDofMapper(m_interfaceStrategy, m_dirStrategy);
    }
    else if ((m_geoTrans == INVERSE_COMPOSITION
              || (m_geoTrans == DIV_CONFORMING
              && m_dirStrategy == dirichlet::nitsche) ) )
    {
        initDofMapper(m_interfaceStrategy, m_dirStrategy);
    }
    else
    {
        GISMO_ERROR("Multipatch with strongly imposed normal DBC for div.pers. splines not implemented!");
    }


    // Initialize the assembler by resize the global matrix and the right-hand side.
    m_matrix.resize(m_idofs,m_idofs);
    m_rhs.resize(m_idofs,1);
    m_matrix.setZero();
    m_rhs.setZero();

    // Estimate max non-zeros per row (only valid for tensor product bases right now)
    // Assumes that degree is equal for all velocity components and patches
    unsigned nzRowsPerColAii = 1;
    unsigned nzRowsPerColB = 1;

    for (int i = 0; i < m_bases[0][0].dim(); ++i)
    {
        nzRowsPerColAii *= 3 * m_bases[0][0].degree(i) + 1;
        nzRowsPerColB   *= 3 * m_bases[tarDim][0].degree(i) + 1; //This is probably not right
    }

    gsVector<index_t> vectorA, vectorB;
    if (m_geoTrans == DIV_CONFORMING)
    {
        vectorA = gsVector<index_t>::Constant(m_idofs - m_dofMapper.back()->freeSize(),
                                          nzRowsPerColAii*tarDim + nzRowsPerColB);
    }
    else
    {
        vectorA = gsVector<index_t>::Constant(m_idofs - m_dofMapper.back()->freeSize(),
                                          nzRowsPerColAii + nzRowsPerColB);
    }
    vectorB = gsVector<index_t>::Constant(m_dofMapper.back()->freeSize(), tarDim*nzRowsPerColAii);
    gsVector<index_t> vectorAB(m_idofs);
    vectorAB << vectorA, vectorB;
    // Reserve space
    // \todo get the fragment of dofs in patch numbered patchIndex
    m_matrix.reserve( vectorAB );
    //gsInfo << "Reserved: " << vectorAB.sum() <<  " and "<< vectorB(0)<< std::endl;

    // Get a block view of the matrix and right hand side:
    //
    // (A11 A12 A13 B1^T)      (rhs_u1)
    // (A21 A22 A23 B2^T)  and (rhs_u2)
    // (A31 A32 A33 B3^T)      (rhs_u3)
    // (B1  B2  B3     0)      (rhs_p )

    //Initialize block structure of m_matrix and m_rhs
    gsVector<index_t> blockPositions(m_dofMapper.size());
    for (unsigned k = 0; k < m_dofMapper.size(); ++k)
    {
        blockPositions[k] = m_dofMapper[k]->freeSize();
    }
    gsVector<index_t> singleCol(1);
    singleCol << 1;

    m_matrixBlocks  = m_matrix.blockView(blockPositions, blockPositions);
    m_rhsBlocks= m_rhs.blockView(blockPositions, singleCol);
}



template<class T>
void gsStokesAssembler<T>::computeDirichletDofs(dirichlet::strategy dirStrategy)
{
    GISMO_ASSERT(!m_bconditions[1], "Pressure Boundary condition not implemented!");

    // Initialize the matrices in m_fixedDofs
    for (unsigned k=0; k< m_fixedDofs.size(); ++k)
    {
        m_fixedDofs[k].resize( m_dofMapper[k]->boundarySize(), 1);
        m_fixedDofs[k].setZero();
    }
    const index_t boundarySize = m_dofMapper[0]->boundarySize();
    const bool interpolate = false;
    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseMatrix<T> globProjMat( boundarySize, boundarySize );
    gsMatrix<T>       globProjRhs;
    globProjRhs.setZero( boundarySize, 1 );

    // Temporaries


    for ( typename gsBoundaryConditions<T>::const_iterator
          it = m_bconditions[0]->dirichletBegin(); it != m_bconditions[0]->dirichletEnd(); ++it )
    {
        //Active_shift is a member of gsGenericBasisEvaluator however I do not have access to it
        unsigned active_shift = 0;

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            continue;
        }

        if (dirStrategy == dirichlet::eliminatNormal)
        {
            index_t norm_basis = (it->ps.side()-1)/2;
            for (int comp = 0; comp < norm_basis; ++comp)
            {
                active_shift += m_bases[comp].size(it->patch());
            }
            if (interpolate)
                computeDirichletDofsComponent(it, norm_basis, active_shift);
            else
                computeDirichletDofsL2Proj(it, norm_basis, globProjMat, globProjRhs, active_shift);

        }
        if (dirStrategy == dirichlet::elimination)
        {
            for (int comp = 0; comp < m_unknownDim[0]; ++comp)
            {
                if (interpolate)
                    computeDirichletDofsComponent(it, comp, active_shift);
                else
                    computeDirichletDofsL2Proj(it, comp, globProjMat, globProjRhs, active_shift);

                active_shift += m_bases[comp].size(it->patch());
            }
        }
    }


    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    if(!interpolate)
    {
        globProjMat.makeCompressed();
        typename gsSparseSolver<T>::CGDiagonal solver;
        m_fixedDofs[0] = solver.compute( globProjMat ).solve ( globProjRhs );
    }
}

template<class T>
void gsStokesAssembler<T>::computeDirichletDofsComponent(
        typename gsBoundaryConditions<T>::const_iterator it,
        index_t comp, unsigned active_shift)
{
    // Get DoFs on this boundary
    gsMatrix<index_t> helper = m_bases[comp][it->patch()].boundary(it->side());
    gsMatrix<index_t> boundary = helper.array() + active_shift;
    //std::cout<< comp<< "Inside computeDirichletDofsComponent number of dofs are " << boundary.rows() << " side FixDdof " << m_fixedDofs[0].rows()<< std::endl;


    // Get the side information
    int dir =  it->side().direction();
    index_t param = ( it->side().parameter() ? 1 : 0);

    // Compute grid of points on the face ("face anchors")
    std::vector< gsVector<T> > rr;
    rr.reserve( m_patches.parDim() );

    for( int i=0; i < m_patches.parDim(); ++i)
    {
        if ( i==dir )
        {
            gsVector<T> b(1);
            b[0] = ( m_bases[comp][it->patch()].component(i).support() ) (0, param);
            rr.push_back(b);
        }
        else
        {
            rr.push_back( m_bases[comp][it->patch()].component(i).anchors().transpose() );
        }
    }

    // Compute Dirichlet values

    //Do invers Piola tarnsform on Dirichlet data
    gsMatrix<T> fpts_row;
    gsMatrix<T> fptsPiolaInv;
    gsMatrix<T> fpts =
        it->function()->eval(m_patches.patch(it->patch()).eval(  gsPointGrid<T>( rr ) ) );
    //if(it->side() == 1 || it->side() == 3)
    //    fpts = fpts;
    if (m_geoTrans == DIV_CONFORMING)
    {
        // Evaluate the geometry
        typename gsGeometryEvaluator<T>::uPtr geoEval(
            getEvaluator(NEED_VALUE|NEED_GRAD_TRANSFORM|NEED_MEASURE, m_patches.patch(it->patch())));
        geoEval->evaluateAt(gsPointGrid<T>( rr ));

        gsMatrix<T> paraPts = gsPointGrid<T>( rr );
        index_t numPts = paraPts.cols();
        fptsPiolaInv.setZero(fpts.rows(), numPts);
        for (index_t kp = 0; kp < numPts; ++kp)
        {
            fptsPiolaInv.col(kp) = (geoEval->gradTransform(kp).transpose())*fpts.col(kp)*geoEval->measure(kp);
        }

        if(it->side() == 1 || it->side() == 3)
            fpts_row = fptsPiolaInv.row(comp);
        else
            fpts_row = fptsPiolaInv.row(comp);

    }
    else
    {
        fpts_row = fpts.row(comp);
    }


    // Interpolate Dirichlet boundary
    typename gsBasis<T>::uPtr h = m_bases[comp][it->patch()].boundaryBasis(it->side());

    typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts_row);
    const gsMatrix<T> & dVals =  geo->coefs();

    // Save corresponding boundary DoFs
    for (index_t j=0; j!= boundary.size(); ++j)
    {
        const int ii= m_dofMapper[0]->bindex( boundary(j) , it->patch() );
        m_fixedDofs[0](ii,0) = dVals(j,0);
    }
    // REMARK: Note that correcting the right-hand-side by the contributions
    // from the fixed DoFs is done/has to be done in the assembling process.
}

// S.Kleiss
/** \brief Incorporates Dirichlet-boundary conditions by L2-projection.
 *
 * ...if the Dirichlet-strategy \em elimination is chosen.\n
 * A global \f$L_2\f$-projection is applied to the given Dirichlet data
 * and the eliminated coefficients are set to the corresponding values.
 * The projection is global in the sense that all Dirichlet-DOFs are
 * computed at once.
 */
template<class T>
void gsStokesAssembler<T>::computeDirichletDofsL2Proj(
        typename gsBoundaryConditions<T>::const_iterator it,
        index_t comp, gsSparseMatrix<T> & globProjMat,
        gsMatrix<T> &globProjRhs, unsigned active_shift)
{
    index_t patchIdx = it->patch();
    const gsBasis<T> & basis = (m_bases[comp])[it->patch()];
    const gsDofMapper & mapper = *m_dofMapper[0];

    typename gsGeometryEvaluator<T>::uPtr geoEval(
        getEvaluator(NEED_VALUE|NEED_GRAD_TRANSFORM|NEED_MEASURE, m_patches.patch(it->patch())));

    // Set up quadrature. The number of integration points in the direction
    // that is NOT along the element has to be set to 1.
    gsVector<index_t> numQuNodes( basis.dim() );
    for( int i=0; i < basis.dim(); i++)
        numQuNodes[i] = (basis.degree(i)+1);
    numQuNodes[  it->side().direction()] = 1;

    gsGaussRule<T> QuRule;
    QuRule.setNodes( numQuNodes);
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    // Create the iterator along the given part boundary.
    typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(it->side());

    for(; bdryIter->good(); bdryIter->next() )
    {
        //bdryIter->computeQuadratureRule( numQuNodes );
        QuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),\
                         quNodes, quWeights);

        geoEval->evaluateAt( quNodes );
        //geoEval->evaluateAt(m_patches[patchIdx].eval( quNodes ));

        gsMatrix<T> fpts = it->function()->eval( m_patches[patchIdx].eval( quNodes ) );

        gsMatrix<T> fpts_row;
        gsMatrix<T> fptsPiolaInv;
        if (m_geoTrans == DIV_CONFORMING)
        {
            gsMatrix<T> paraPts = m_patches[patchIdx].eval( quNodes );
            index_t numPts = paraPts.cols();
            fptsPiolaInv.setZero(fpts.rows(), numPts);
            for (index_t kp = 0; kp < numPts; ++kp)
            {
                fptsPiolaInv.col(kp) = (geoEval->gradTransform(kp).transpose())*fpts.col(kp)*geoEval->measure(kp);
            }

            fpts_row = fptsPiolaInv.row(comp);

        }
        else
        {
            fpts_row = fpts.row(comp);
        }

        gsMatrix<T> basisVals;
        basis.eval_into( quNodes, basisVals);

        // Indices involved here:
        // --- Local index:
        // Index of the basis function/DOF on the patch.
        // Does not take into account any boundary or interface conditions.
        // --- Global Index:
        // Each DOF has a unique global index that runs over all patches.
        // This global index includes a re-ordering such that all eliminated
        // DOFs come at the end.
        // The global index also takes care of glued interface, i.e., corresponding
        // DOFs on different patches will have the same global index, if they are
        // glued together.
        // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
        // The eliminated DOFs, which come last in the global indexing,
        // have their own numbering starting from zero.

        // Get the global indices (second line) of the local
        // active basis (first line) functions/DOFs:
        gsMatrix<index_t> globIdxAct;
        basis.active_into(quNodes.col(0), globIdxAct );
        globIdxAct = globIdxAct.array() + active_shift;
        mapper.localToGlobal( globIdxAct, patchIdx, globIdxAct);

        // Out of the active functions/DOFs on this element, collect all those
        // which correspond to a boundary DOF.
        // This is checked by calling mapper::is_boundary_index( global Index )

        // eltBdryFcts stores the row in basisVals/locIdxAct, i.e.,
        // something like a "element-wise index"
        std::vector<index_t> eltBdryFcts;
        for( index_t i=0; i < globIdxAct.rows(); i++)
            if( mapper.is_boundary_index( globIdxAct(i,0)) )
                eltBdryFcts.push_back( i );

        // Do the actual assembly:
        for( index_t k=0; k < quNodes.cols(); k++ )
        {
            T weight_k = quWeights[k] * geoEval->measure(k);;

            // Only run through the active boundary functions on the element:
            for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
            {
                // Each active boundary function/DOF in eltBdryFcts has...
                // ...the above-mentioned "element-wise index"
                unsigned i = eltBdryFcts[i0];
                // ...the boundary index.
                unsigned ii = mapper.global_to_bindex( globIdxAct( i ));
                for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                {
                    unsigned j = eltBdryFcts[j0];
                    unsigned jj = mapper.global_to_bindex( globIdxAct( j ));

                    // Use the "element-wise index" to get the needed
                    // function value.
                    // Use the boundary index to put the value in the proper
                    // place in the global projection matrix.
                    globProjMat.coeffRef(ii,jj) += weight_k * basisVals(i,k) * basisVals(j,k);
                } // for j

                globProjRhs(ii,0) += weight_k *  (basisVals(i,k) * fpts_row.col(k).transpose()).value();;

            } // for i
        } // for k
    } // bdryIter

} // computeDirichletDofsL2Proj

// Assembles the final system with all boundary conditions contained
template<class T>
void gsStokesAssembler<T>::assemble()
{

    // If the Dirichlet strategy is elimination then precompute
    // Dirichlet dofs (m_dofMapper excludes these from the system)
    if ( m_dirStrategy == dirichlet::elimination || 
         m_dirStrategy == dirichlet::eliminatNormal)
    {
        computeDirichletDofs(m_dirStrategy);
    }


    if (m_idofs == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed dirichlet boundary only.\n" <<"\n" ;
        return;
    }
    // Assemble the system matrix and right-hand side
    else
    {
        for (size_t patchIndex=0; patchIndex < m_patches.nPatches(); ++patchIndex )
        {
            // Assemble stiffness matrix and rhs for the
            // local patch and add to m_matrix and m_rhs
            assemblePatch(patchIndex);
        }
    }

    // Enforce Dirichlet boundary conditions by Nitsche's method
    if ( m_dirStrategy == dirichlet::nitsche ||
         m_dirStrategy == dirichlet::eliminatNormal)//
    {
        for ( typename gsBoundaryConditions<T>::const_iterator
              BCiterator = m_bconditions[0]->dirichletBegin();
              BCiterator != m_bconditions[0]->dirichletEnd(); ++BCiterator )
        {
            gsStokesAssembler<T>::boundaryNitsche(BCiterator->patch(),
                                               *BCiterator->function(), BCiterator->side());
            //boundaryNitschePressure(BCiterator->patch(), *BCiterator->function(), BCiterator->side());
            //freeAll( basis_vector);
        }
    }

    // Enforce Neumann boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          BCiterator = m_bconditions[0]->neumannBegin();
          BCiterator != m_bconditions[0]->neumannEnd(); ++BCiterator )
    {
        gsStokesAssembler<T>::boundaryNeumann(BCiterator->patch(),
                                           *BCiterator->function(),BCiterator->side() ); 
        //freeAll( B_vec);
    }


    // Enforce Pressure boundary conditions
    if (m_bconditions[1]) // If pressure BC
    {
        GISMO_ERROR("Dirichlet condition for pressure not implemented as it alone is not sufficient to get coersivity. Use instead grad(u)*n + pn = h(x) as a Neumann condition on u.");
        // Changeset 2820 as some code trying to add Dirichlet pressure BC.
    }
    m_matrix.makeCompressed();

}


template<class T>
void gsStokesAssembler<T>::assemblePatch(int patchIndex)
{
    //Comment explinations:
    //phi_j1 are the basis functions for the velocity field's first component
    //phi_j2 are the basis functions for the velocity field's second component
    //phi_j3 are the basis functions for the velocity field's thired component
    //psi_j  are the basis functions for the pressue

    //Target dimention for the velocity
    const short_t tarDim = m_unknownDim[0];
    const int tar2      = tarDim*tarDim;

    const T nu = m_nu;

    // Copy data
    std::vector<gsBasis<T> *> basis_u_vec;
    for (int k = 0; k < tarDim; ++k)
    {
        basis_u_vec.push_back(&m_bases[k][patchIndex]);
    }
    gsBasis<T> & basis_p = m_bases.back()[patchIndex];

    // Quadrature  nodes
    gsVector<index_t> numNodes(tarDim);
    for (index_t comp = 0 ; comp < tarDim; ++comp)
    {
        numNodes(comp) = basis_u_vec[0]->degree(comp) + 1;
    }

    gsGaussRule<T> quad(numNodes);
    gsMatrix<T> quNodes  ; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights

    // values of the right-hand side
    gsMatrix<T> rhsVals;
    gsMatrix<T> localMatrixA;
    gsMatrix<T> localMatrixB;

    // position of top-left corner of B matrix
    const index_t tShift = 0;
    const index_t pShift = m_dofMapper[0]->freeSize();

    // local load vector
    gsVector<T> localRhs_p;
    gsMatrix<T> localRhs_u;

    // Evaluate the geometry
    typename gsGeometryEvaluator<T>::uPtr geoEval(
        getEvaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER, m_patches.patch(patchIndex)));


    // Make domain element iterator
    typename gsBasis<T>::domainIter domIter = basis_u_vec[0]->makeDomainIterator();

    //Initilize basis evaluator
    gsBasisEvaluator<T>     &eval_u=*makeBasisEvaluator( basis_u_vec,
                                                               NEED_VALUE | NEED_GRAD | NEED_DIV,
                                                               &m_patches.patch(patchIndex),
                                                               m_geoTrans);
    gsBasisEvaluator<T>     &eval_p=*makeBasisEvaluator( basis_p,
                                                               NEED_VALUE | NEED_GRAD,
                                                               &m_patches.patch(patchIndex),
                                                               INVERSE_COMPOSITION);

    gsSparseMatrix<T> rhs_mod;

    //Initialize the local to global methods
    gsL2GMapped<>            vel_mappers(m_matrix, rhs_mod, *m_dofMapper[0], *m_dofMapper[0], patchIndex);
    gsL2GMappedMultiplier<>  div_mappers(m_matrix, rhs_mod, *m_dofMapper[0], *m_dofMapper[1], patchIndex, tShift, pShift);
    gsL2GMappedRhs<>         rhs_mappers(m_rhs, *m_dofMapper[0], patchIndex);


    // this would more elegantly be in another loop
    rhs_mod.resize(m_matrix.rows(),m_fixedDofs[0].rows());

    // Start iteration over elements
    for (; domIter->good(); domIter->next())
    {
        // Map the Quadrature rule to the element
        quad.mapTo( domIter->lowerCorner(),
                    domIter->upperCorner(), quNodes, quWeights);

        // Evaluate basis functions on element
        geoEval->evaluateAt(quNodes);

        eval_u.evaluateAt (quNodes, *geoEval);
        eval_p.evaluateAt (quNodes, *geoEval);

        const gsMatrix<T> &basisValues_u(eval_u.values()); //Transformed basis values
        const gsMatrix<T> &basisGrad_u(eval_u.derivs());   //Transformed gardients values
        const gsMatrix<T> &basisDiv_u(eval_u.divs());      //Transformed divergence values
        const gsMatrix<T> &basisValue_p(eval_p.values());  //Transformed values

        // Active basis functions at one quadrature node;
        const gsMatrix<index_t> &actives_u(eval_u.actives());
        const gsMatrix<index_t> &actives_p(eval_p.actives());

        index_t numActU = actives_u.rows(); //Number of Active components for all velocity components
        index_t numActp = actives_p.rows(); //Number of Active components for pressure

        // Evaluate right-hand side at the geometry points
         m_rhs_function->eval_into( geoEval->values(), rhsVals );


        localMatrixA.setZero(numActU, numActU);//local_A
        localMatrixB.setZero(numActp, numActU);//local_B
        localRhs_u.setZero(numActU, tarDim);
        localRhs_p.setZero(numActp);


        // Loop over quadrature nodes for velocity
        for (index_t node = 0; node < quNodes.cols(); ++node)
        {
            // weight * abs(det J), where J is geometry Jacobian.
            T weight = quWeights(node) * geoEval->measure(node);

            for (index_t ai = 0; ai< numActU; ai++)
            {
                for (index_t aj = 0; aj< numActU; aj++)
                {
                    //Siffness matrix for A_ij: nu (D phi, D phi)
                    //Frobenius mutiplication
                    localMatrixA(ai, aj) += weight * nu *
                            (basisGrad_u.block(tar2*(ai) ,node, tar2,1).transpose() *
                             basisGrad_u.block(tar2*(aj) ,node, tar2,1) ).value();
                }
                // B_i: (div(phi), psi)
                for (index_t ap = 0; ap< numActp; ap++)
                {
                    localMatrixB(ap,ai) +=
                            weight * basisDiv_u(ai,node) * basisValue_p(ap,node);
                }
                //Calculateing the b_i local load vector: (phi,f)
                localRhs_u(ai,0) += weight *
                        (basisValues_u.block(tarDim*(ai) ,node, tarDim,1).transpose() *
                         rhsVals.col(node)).value();
            }
        }  // end loop Gauss nodes

        //----- ADD LOCAL TO GLOBAL -----//
        vel_mappers.store(actives_u, actives_u, localMatrixA);
        div_mappers.store(actives_u, actives_p, localMatrixB.transpose());
        rhs_mappers.store(actives_u, actives_u, localRhs_u.col(0));

    } // end loop over all domain elements


    // add rhs modifications
    m_rhs-=rhs_mod*(m_fixedDofs[0].col(0));

    delete &eval_u;
    delete &eval_p;
}

template<class T>  void
gsStokesAssembler<T>::boundaryNeumann( const int patchIndex,
                                       const gsFunction<T> & f,
                                       const boxSide s)
{
    // Copy data
    const index_t tarDim = m_unknownDim[0];
    std::vector<gsBasis<T> *> basis_u_vec;
    for (int k = 0; k < tarDim; ++k)
    {
        basis_u_vec.push_back(&m_bases[k][patchIndex]);
    }

    // Quadrature  nodes
    gsVector<index_t> numNodes(tarDim);
    for (index_t comp = 0 ; comp < tarDim; ++comp)
    {
        if (comp == s.direction())
            numNodes(comp) =  1;
        else
            numNodes(comp) = basis_u_vec[0]->degree(comp) + 1;

    }
    gsGaussRule<T> quad(numNodes);

    gsMatrix<T> quNodes  ; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights

    typename gsGeometryEvaluator<T>::uPtr
            geoEval(getEvaluator(NEED_VALUE | NEED_MEASURE, m_patches.patch(patchIndex)));

    // Temporaries
    gsMatrix<T> fev, basisValues;
    gsVector<T> unormal(tarDim);


    // Active basis functions at one quadrature node
    gsMatrix<index_t> actives_u;

    // Local load vactor
    gsMatrix<T> localRhs_u;

    // Make domaboxSiderator
    typename gsBasis<T>::domainIter domIter = basis_u_vec[0]->makeDomainIterator(s);

    //Initilize basis evaluator
    gsBasisEvaluator<T>     &basis_eval_u=*makeBasisEvaluator(
         basis_u_vec,NEED_VALUE, &m_patches.patch(patchIndex), m_geoTrans);

    typedef gsShiftWriter<gsMatrix<T> >      SWR;
    SWR shift_rhs(m_rhs, 0, 0);
    gsLocalToGlobalMapper<T> * rhs_mappers = new gsL2GMappedRhs<SWR> (shift_rhs, *m_dofMapper[0], patchIndex);

    // Iterate over all boundary grid cells
    for (; domIter->good(); domIter->next())
    {
        // Map the Quadrature rule to the element
        quad.mapTo( domIter->lowerCorner(),
                    domIter->upperCorner(), quNodes, quWeights);


        // Evaluate the geometry and basis functions on element
        geoEval->evaluateAt(quNodes);
        basis_eval_u.evaluateAt (quNodes, *geoEval);

        actives_u = basis_eval_u.actives();

        //ev and basisGrads
        basisValues = basis_eval_u.values(); //Transformed velocity values

        index_t numActU = actives_u.rows(); //Number of Active components for all velocity components

        // Evaluate the Dirichlet data
        f.eval_into(geoEval->values(), fev);

        // Initialize local rhs to 0
        localRhs_u.setZero(numActU, 1);


        // Loop over quadrature nodes
        for (index_t node = 0; node < quNodes.cols(); ++node)
        {
            geoEval->outerNormal(node, s, unormal);

            T weight = quWeights(node) * unormal.norm();

            //Calculateing the right hand side
            for (index_t ai = 0; ai< numActU; ai++)
            {
                localRhs_u(ai,0) += weight *
                        (fev.col(node).transpose() *
                         basisValues.block(tarDim*ai ,node, tarDim,1)).value();
            }
        }

        rhs_mappers->store(actives_u, actives_u, localRhs_u.col(0));

    } // end loop over all domain elements

    // clean up other stuff
    delete rhs_mappers;
    delete &basis_eval_u;
}

template<class T>  void
gsStokesAssembler<T>::boundaryNitsche( const int patchIndex,
                                       const gsFunction<T> & f,
                                       const boxSide s)
{
    // Copy data
    const index_t tarDim = m_unknownDim[0];
    const index_t tar2 = tarDim*tarDim;
    std::vector<gsBasis<T> *> basis_u_vec;
    for (int k = 0; k < tarDim; ++k)
    {
        basis_u_vec.push_back(&m_bases[k][patchIndex]);
    }
    const gsBasis<T> & basis_p = m_bases.back()[patchIndex];
    const T kappa = m_nu;

    //const T mu = gsAssemblerUtils<T>::getMu(*basis_u_vec[0]);
    T mu = 5*(basis_u_vec[1][0].degree(0) + 1);//Evans


    //const std::vector<gsDofMapper  *> & mapper = m_dofMapper;

    // Quadrature  nodes
    gsVector<index_t> numNodes(tarDim);
    for (index_t comp = 0 ; comp < tarDim; ++comp)
    {
        if (comp == s.direction())
            numNodes(comp) =  1;
        else
            numNodes(comp) = basis_u_vec[0]->degree(comp) + 1;

    }
    gsGaussRule<T> quad(numNodes);
    gsMatrix<T> quNodes  ; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights

    unsigned flags;
    if (m_geoTrans == DIV_CONFORMING)
        flags = NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM | NEED_2ND_DER | NEED_MEASURE;
    else
        flags = NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM | NEED_MEASURE;
    typename gsGeometryEvaluator<T>::uPtr
            geoEval(getEvaluator(flags, m_patches.patch(patchIndex)));

    // Temporaries
    gsMatrix<T> fev, tmp_gradTj, tmp_gradTi;
    gsVector<T> unormal(tarDim), tmp_gradNormi, tmp_gradNormj;

    // basis(Values/Grad) contains the stacked (Values/Gradients)
    // of all basis functions at one quadrature node for each collum
    gsMatrix<T> basisValues, basisGrad, pressureValues;

    // Active basis functions at one quadrature node
    gsMatrix<index_t> actives_u;
    gsMatrix<index_t> actives_p;


    // (dense) local stiffness matrice within one grid cell
    gsMatrix<T> localMatrixA;
    gsMatrix<T> localMatrixB;

    // Local load vactor
    gsMatrix<T> localRhs_u;
    gsVector<T> localRhs_p;

    //size of momentum matrix(A)
    const index_t pShift = m_dofMapper[0]->freeSize();;

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIter = basis_u_vec[(s-1)/2]->makeDomainIterator(s);



    //Initilize basis evaluator
    gsBasisEvaluator<T>     &basis_eval_u=*makeBasisEvaluator( basis_u_vec,
                                                               NEED_VALUE | NEED_GRAD | NEED_DIV,
                                                               &m_patches.patch(patchIndex),
                                                               m_geoTrans);
    gsBasisEvaluator<T>     &basis_eval_p=*makeBasisEvaluator( basis_p,
                                                               NEED_VALUE | NEED_GRAD,
                                                               &m_patches.patch(patchIndex),
                                                               INVERSE_COMPOSITION);

    //Initialize the local to global methods
    gsLocalToGlobalMapper<T> *vel_mappers = new gsLocalToGlobalMapper<T>;
    gsLocalToGlobalMapper<T> *div_mappers = new gsLocalToGlobalMapper<T>;
    gsLocalToGlobalMapper<T> *rhs_mappers_u = new gsLocalToGlobalMapper<T>;
    gsLocalToGlobalMapper<T> *rhs_mappers_p = new gsLocalToGlobalMapper<T>;

    gsSparseMatrix<T> rhs_mod;
    gsSparseMatrix<T> rhs_mod1;


    typedef gsShiftWriter<gsSparseMatrix<> > SW;
    typedef gsShiftWriter<gsMatrix<T> >      SWR;

    unsigned tShift=0;
    unsigned uShift=0;

    SW shift_m(m_matrix, tShift, uShift);
    SW shift_rhsMod(rhs_mod, tShift, 0);
    
    delete (vel_mappers); // release memory before assigning a new value
    vel_mappers= new gsL2GMapped<SW,SW>(
                shift_m, shift_rhsMod,
                *m_dofMapper[0], *m_dofMapper[0],
                patchIndex);

    // this would more elegantly be in another loop
    rhs_mod.resize(m_matrix.rows(),m_fixedDofs[0].rows());
    rhs_mod1.resize(m_matrix.rows(),m_fixedDofs[0].rows());

    delete (div_mappers); // release memory before assigning a new value
    div_mappers = new gsL2GMappedMultiplier<gsSparseMatrix<T> >(
                m_matrix, rhs_mod1,
                *m_dofMapper[0], *m_dofMapper[1],
                patchIndex, tShift, pShift);

    gsMatrix<T>       m_rhs_test;
    m_rhs_test.setZero(m_rhs.rows(), m_rhs.cols());

    SWR shift_rhs_u(m_rhs, tShift, 0);
    delete (rhs_mappers_u); // release memory before assigning a new value
    rhs_mappers_u = new gsL2GMappedRhs<SWR> (shift_rhs_u, *m_dofMapper[0], patchIndex);
    SWR shift_rhs_p(m_rhs_test, pShift, 0);
    delete (rhs_mappers_p); // release memory before assigning a new value
    rhs_mappers_p = new gsL2GMappedRhs<SWR> (shift_rhs_p, *m_dofMapper[1], patchIndex);


    // Iterate over all boundary grid cells
    for (; domIter->good(); domIter->next())
    {
        // Map the Quadrature rule to the element
        quad.mapTo( domIter->lowerCorner(),
                                      domIter->upperCorner(), quNodes, quWeights);


        // Evaluate the geometry and basis functions on element
        geoEval->evaluateAt(quNodes);

        basis_eval_u.evaluateAt (quNodes, *geoEval);
        basis_eval_p.evaluateAt (quNodes, *geoEval);

        actives_u = basis_eval_u.actives();
        actives_p = basis_eval_p.actives();


        //ev and basisGrads
        basisValues = basis_eval_u.values(); //Transformed velocity values
        basisGrad   = basis_eval_u.derivs(); //Transformed gardients values
        pressureValues = basis_eval_p.values();

        index_t numActU = actives_u.rows(); //Number of Active components for all velocity components
        index_t numActp = actives_p.rows(); //Number of Active components for pressure


        //Find cell size.
        T h_Q = domIter->getPerpendicularCellSize();
        T jac_inf = geoEval->jacobians().maxCoeff();
        if ( math::abs(geoEval->jacobians().minCoeff()) > jac_inf )
            jac_inf = math::abs(geoEval->jacobians().minCoeff());
        jac_inf = 1.0;
        T h_K = h_Q*jac_inf;
        //gsDebug<< "jac_inf: " << jac_inf << "  h_Q: "<< h_Q<< "  h_K: "<< h_K<< "  nu: "<< mu << "  mu: "<< mu/h_K<< " side: "<< s <<"\n";

        // Evaluate the Dirichlet data
        f.eval_into(geoEval->values(), fev);

        // Initialize local linear system to 0
        localMatrixA.setZero(numActU, numActU);//local_A
        localMatrixB.setZero(numActp, numActU);//local_B
        localRhs_u.setZero(numActU, tarDim);
        localRhs_p.setZero(numActp);


        // Loop over quadrature nodes
        for (index_t node = 0; node < quNodes.cols(); ++node)
        {
            //gsDebug <<"node: "<<node<<"\n";
            // Compute the outer normal vector
            geoEval->outerNormal(node, s, unormal);

            T weight = kappa * quWeights(node) * unormal.norm();

            // Compute the unit normal vector
            unormal.normalize();
            //Calculateing the right hand side
            for (index_t ai = 0; ai< numActU; ai++)
            {
                tmp_gradTi.noalias() = basisGrad.block(tar2*(ai) ,node, tar2,1);
                tmp_gradTi.resize(tarDim,tarDim);
                tmp_gradNormi.noalias() = unormal.transpose()*tmp_gradTi;
                localRhs_u(ai,0) += weight *
                        ((mu/h_K)*(fev.col(node).transpose() *
                                   basisValues.block(tarDim*(ai) ,node, tarDim,1)).value()
                         - fev.col(node).dot((tmp_gradNormi)));
                // A: (grad(u)*n, v)
                for (index_t aj = 0; aj< numActU; aj++)
                {   //JS2: This can be optimised because of symmertri properties!
                    tmp_gradTj = basisGrad.block(tar2*(aj) ,node, tar2,1);
                    tmp_gradTj.resize(tarDim,tarDim);

                    tmp_gradNormj.noalias() = unormal.transpose()*tmp_gradTj;
                    localMatrixA(ai,aj) += weight *
                            ((mu/h_K)*(basisValues.block(tarDim*(aj) ,node, tarDim,1).transpose() *
                                       basisValues.block(tarDim*(ai) ,node, tarDim,1)).value() -
                             (basisValues.block(tarDim*(ai) ,node, tarDim,1).transpose() * tmp_gradNormj).value() -
                             (basisValues.block(tarDim*(aj) ,node, tarDim,1).transpose() * tmp_gradNormi).value());
                }
                // B: (div(phi), psi)
                for (index_t ap = 0; ap< numActp; ap++)
                {
                    localMatrixB(ap,ai) -= weight/kappa * pressureValues(ap,node) *
                            (unormal.transpose() * basisValues.block(tarDim*(ai) ,node, tarDim,1)).value();
                }
                T testtmp = (unormal.transpose() * basisValues.block(tarDim*(ai) ,node, tarDim,1)).value();
                if (testtmp<0){testtmp *= -1;}
                if(m_dofMapper[0]->is_free_index(m_dofMapper[0]->index( actives_u(ai,0) , 0 )))
                {
                    if (testtmp > 1e-12)
                    {
                        //std::cout << "Local " << ai << " global: " << actives_u(ai,0) << " GLOBAL: " << m_dofMapper[0]->is_free_index(m_dofMapper[0]->index( actives_u(ai,0) , 0 )) << " side " << s <<" normal: " << unormal.transpose()<< " values: " << basisValues.block(tarDim*(ai) ,node, tarDim,1).transpose() << std::endl;
                        //if(! m_dofMapper[0]->is_free_index(m_dofMapper[0]->index( actives_u(ai,0) , 0 )))
                        //{
                        //    std::cout << "Jarle " << std::endl;
                        //}
                    }
                    else
                        {}//std::cout << "Under TOT" << std::endl;
                }
                else
                {}//std::cout << "Local " << ai << " global: " << actives_u(ai,0) <<  "Global " <<  m_dofMapper[0]->bindex( actives_u(ai,0) , 0 ) <<" side " << s <<" BOUNDARY!!!" << std::endl;
            }
            for (index_t ap = 0; ap< numActp; ap++)
            {
                localRhs_p(ap) -= weight/kappa * unormal.dot(fev.col(node)) * pressureValues(ap,node);
            }
        }// end loop Gauss nodes

        //----- ADD LOCAL TO GLOBAL -----//
        vel_mappers->store(actives_u, actives_u, localMatrixA);
        rhs_mappers_u->store(actives_u, actives_u,localRhs_u.col(0));
        //If the normal component of the Dirichlet BC is strongly forced,
        //don't add pressure contribution.
        //NB! If projection is used to find boundary DOFs,
        //than these contribution must be added to get correct accuracy!
        if (!(m_dirStrategy == dirichlet::eliminatNormal))
        {
            div_mappers->store(actives_u, actives_p, localMatrixB.transpose()); //rhs_nod1
            rhs_mappers_p->store(actives_p, actives_p,localRhs_p); //m_rhs_test
        }

    } // end loop over all domain elements

    // add rhs modifications
    m_rhs-=rhs_mod*(m_fixedDofs[0].col(0));
    //NB! If projection is used to find boundary DOFs,
    //than these contribution must be added to get correct accuracy!
    if (!(m_dirStrategy == dirichlet::eliminatNormal))
    {
        m_rhs -= rhs_mod1*(m_fixedDofs[0].col(0)); //(p*n,u)_\Gamma contributions
        m_rhs += m_rhs_test;
    }
    //std::cout <<s<<"  SIDE   " << std::endl;

    for (index_t k = 0;k <m_rhs_test.rows();++k)
    {
        int tmp_s = s;

        if (tmp_s == 1 || true)
        {
            //T atmp = m_rhs_test(k,0);
            //T btmp = rhs_mod1.row(k)*(m_fixedDofs[0].col(0));
            //if (std::abs(atmp-btmp)> 3e-15)
            //    std::cout <<k<<"=  k   " << atmp << "   "  << btmp<< " differnece "<<atmp-btmp<< std::endl;
        }
    }
    //std::cout << rhs_mod1*(m_fixedDofs[0].col(0))<< std::endl;

    // free mappers
    delete vel_mappers;
    delete div_mappers;
    delete rhs_mappers_u;
    delete rhs_mappers_p;

    // clean up other stuff
    delete &basis_eval_u;
    delete &basis_eval_p;
}



// Solves the linear system and fills up \a m_sysSolution
template<class T>
void gsStokesAssembler<T>::solveSystem()
{

    std::cout << "Solve linear system of size: " << m_idofs << "\n";

    //Direct solvers
    //typename gsSparseSolver<T>::LU solver;
    //typename gsSparseSolver<T>::BiCGSTABILUT solver;
    typename gsSparseSolver<T>::QR solver; //COLAMDOrdering<int>
    gsInfo << "Using direct solver" << std::endl;

    //solver.compute( m_matrix );
    solver.analyzePattern(m_matrix);
    solver.factorize(m_matrix);
    gsDebug << "   Eigen  lastErrorMessage: " << solver.lastErrorMessage () << "\n";
    m_sysSolution = solver.solve(m_rhs);
    //_sysSolution = solver.compute(m_matrix).solve (m_rhs);

    gsInfo << "Solved with SparceLU\n" ;
    //gsInfo << "Solved with BiCGSTAB \n" ;
}

// Solution field(s)
template<class T>
gsFunction<T> * gsStokesAssembler<T>::reconstructPatchSolution(index_t unk, int p, bool hasMatrixRhs) const
{
    gsDebug << "Im now using the derived class reconstruction\n";
    const gsMatrix<T> & data = m_sysSolution;

    GISMO_ASSERT ( (signed) m_bases.size() == m_unknownDim.sum(), "Implementation assumes one basis for each unknown component");

    // Target dimension for unknown unk
    const index_t tarDim = m_unknownDim[unk];
    // The position of the unknown in the m_sysSolution vector/matrix
    int UnkPosSysSol = 0;

    //Find the basis index
    index_t basis_ind = 0;
    for (index_t k = 0 ; k < unk; ++k)
        basis_ind += m_unknownDim[k];

    std::vector<index_t> patchShifts(tarDim);
    for (index_t comp = 0; comp < tarDim; ++comp)
    {
        patchShifts[comp] = 0;
        for (index_t k = 0; k < p; ++k)
        {
        patchShifts[comp] += m_bases[basis_ind + comp][k].size();
        }

    }
    std::vector<index_t> basesShifts(tarDim);
    basesShifts[0] = 0;
    for (index_t k = 1; k < tarDim; ++k)
    {
        basesShifts[k] = basesShifts[k-1] + m_bases[basis_ind + k - 1].totalSize();
    }

    if (unk==0)
    {
        gsDebug << "size of Basis 0: " << m_bases[0][0].size()<< "\n";
        gsDebug << "Degree of Basis 0: " << m_bases[0][0].degree(0)<< " and "<< m_bases[0][0].degree(1)<<"\n";
        gsDebug << "size of Basis 1: " << m_bases[1][0].size()<< "\n";
        gsDebug << "Degree of Basis 1: " << m_bases[1][0].degree(0)<< " and "<< m_bases[1][0].degree(1)<<"\n";
        gsDebug << "size of Basis 2: " << m_bases[2][0].size()<< "\n";
        gsDebug << "Degree of Basis 2: " << m_bases[2][0].degree(0)<< " and "<< m_bases[2][0].degree(1)<<"\n";
    }

    for (int k = 0; k < unk; ++k)
    {
        // Finding position in solution vector (for unknowns indexes less then unk)
        UnkPosSysSol += m_dofMapper[k]->freeSize();
    }
    index_t szCoeff = 0;
    for (index_t k = 0; k < tarDim; ++k)
    {
        szCoeff += m_bases[basis_ind + k][p].size();
    }
    if (unk == 0 && m_geoTrans == DIV_CONFORMING)
    {
        gsMatrix<T> coeffs(szCoeff, 1);
        index_t cShift = 0;
        for (index_t comp = 0; comp < tarDim; ++comp)
        {
            for (index_t i = 0; i < m_bases[basis_ind + comp][p].size(); ++i)
            {
                if ( m_dofMapper[unk]->is_free(i + basesShifts[comp], p) ) // internal or interface
                {
                    coeffs(i + cShift, 0) = data(UnkPosSysSol + m_dofMapper[unk]->index(i + basesShifts[comp], p),0);
                }
                else // eliminated DoFs: fill with Dirichlet data
                {
                    coeffs(i + cShift, 0) = m_fixedDofs[unk]( m_dofMapper[unk]->bindex(i + basesShifts[comp], p),0);
                }
                // Increase position with the size of the component
            }
            cShift += m_bases[basis_ind + comp][p].size();
        }
        std::vector<gsBasis<T> *> bases_vec;
        for (index_t k = 0; k < tarDim; ++k)
        {
            bases_vec.push_back( const_cast<gsBasis<T> *>(&m_bases[basis_ind + k][p]) );
        }
        return new gsDivConSolution<T> (coeffs,m_patches[p],bases_vec);

    }
    else
    {
        //Coefficient matrix
        gsMatrix<T> coeffs(m_bases[basis_ind][p].size(), tarDim);

        for (index_t comp = 0; comp < tarDim; ++comp)
        {
            for (index_t i = 0; i < m_bases[basis_ind + comp][p].size(); ++i)
            {
                if ( m_dofMapper[unk]->is_free(i + basesShifts[comp], p) ) // internal or interface
                {
                    coeffs(i,comp) = data(UnkPosSysSol + m_dofMapper[unk]->index(i + basesShifts[comp], p),0);
                }
                else // eliminated DoFs: fill with Dirichlet data
                {
                    coeffs(i,comp) = m_fixedDofs[unk]( m_dofMapper[unk]->bindex(i + basesShifts[comp], p),0);
                }
                // Increase position with the size of the component
            }
        }

        //If reconstructing the pressure:
        if (unk == 1)
            return m_bases[basis_ind][p].makeGeometry( give(coeffs) ).release();
        //If reconstruction the velocity:
        else if (unk == 0)
            return m_bases[0][p].makeGeometry( give(coeffs) ).release();
        else
            GISMO_ERROR("Number of unknowns is wrong");
    }
}

} // namespace gismo


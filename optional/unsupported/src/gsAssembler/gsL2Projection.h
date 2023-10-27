/** @file gsL2Projection.h

    @brief Computes the L2 projection.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Beiser
*/


#pragma once

namespace gismo
{

/** @brief The gsL2Projection class provides the functionality
 * to calculate the L2 Projection of a gsFunction to a gsField.
 *
 * \ingroup Assembler
*/

template <class T>
class gsL2Projection
{
// Constructor
public:

    gsL2Projection(const gsFunction<T> & _Function,
             const gsMultiPatch<T> & _MultiPatch,
             const gsMultiBasis<T> & _MultiBasis,
             gsMultiPatch<> & _ProjectionMultiPatch) 
    :   Function( &_Function )
    ,   MultiPatch( &_MultiPatch )
    ,   MultiBasis( &_MultiBasis )
    ,   ProjectionMultiPatch( &_ProjectionMultiPatch )
    {

        // Initialize and define assembler for L2 projection 
        gsExprAssembler<> ProjectionExprAssembler(1,1);
        ProjectionExprAssembler.setIntegrationElements( *MultiBasis );
        gsExprAssembler<>::geometryMap ProjectionG = ProjectionExprAssembler.getMap( *MultiPatch );
        gsExprAssembler<>::space       ProjectionH = ProjectionExprAssembler.getSpace( *MultiBasis );
        gsExprAssembler<>::variable    Projectionf = ProjectionExprAssembler.getCoeff( *Function, ProjectionG );

        // Assembler mass matrix and rhs 
        ProjectionExprAssembler.initSystem();
        ProjectionExprAssembler.assemble( ProjectionH * ProjectionH.tr() * meas(ProjectionG), ProjectionH * Projectionf * meas(ProjectionG) );

        // Initialize solver and solve system
        gsSparseSolver<>::BiCGSTABDiagonal SparseSolver;
        SparseSolver.compute( ProjectionExprAssembler.matrix() );
        gsMatrix<> Projectionx = SparseSolver.solve( ProjectionExprAssembler.rhs() );

        // Initialize solution and transfer to gsMultiPatch
        gsExprAssembler<>::solution ProjectionSolution = ProjectionExprAssembler.getSolution( ProjectionH, Projectionx );
        ProjectionSolution.extract( *ProjectionMultiPatch );

    }

// member functions
public:
    
    gsField<T> computeField()
    {
        gsField<> ProjectionField( *MultiPatch, *ProjectionMultiPatch );
        return ProjectionField;
    }

    gsMultiPatch<T> computeMultiPatch()
    {
        return *ProjectionMultiPatch;
    }

// member variables
private: 

    const gsFunction<T>* Function;
    const gsMultiPatch<T>* MultiPatch;
    const gsMultiBasis<T>* MultiBasis;
    gsMultiPatch<T>* ProjectionMultiPatch;

};


} // namespace gismo
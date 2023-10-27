/**  gsPoissonCookBook.h

    Ingredients for the poisson equation for the recipe assembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
*/

#pragma once

#include <gsRecipeAssembler/gsRecipe.h>
#include <gsRecipeAssembler/gsPDEOperators.h>
#include <gsRecipeAssembler/gsLocalToGlobal.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsMath.h>

#include<vector>

/** \file here we collect some utils to solve Poisson using the gsRecipeAssemble function
    we assume that the unknown space and the test space are the same, but the general case
    is obtained by changing the test space field (tSpace) of the RecipeEntries

    the other limitation is that the unknown function must be scalar valued.
    this is because for scalar valued functions we can provide a vector valued rhs functions
    meaning we want to solve different problems.
    For vector valued functions this is not possible and a different operator is needed for the test.
**/

namespace gismo {


struct gsPoissonCookBook {
    // ingredients for volume integrals
    enum
    {
        STIFFNESS = 0,
        RHS = 1
    };
    template <typename T> // to do: rename to getRecipe
    static gsRecipe<T> makePoissonRecipe (const gsFunction<T>  &rhsD, gsSparseMatrix<T> &stiff, gsMatrix<T> &rhs)
    {
        gsRecipe<T> result(2);
        result[STIFFNESS].setOperator(new gsGradGradOp<T>());
        result[STIFFNESS].setTestSpace(0);
        result[STIFFNESS].setUnknownSpace(0);
        result[STIFFNESS].setRule(new gsL2GActives<gsSparseMatrix<T> >(stiff));

        result[RHS].setOperator(new gsL2TestOp<T>(rhsD));
        result[RHS].setTestSpace(0);
        result[RHS].setUnknownSpace(0);
        result[RHS].setRule(new gsL2GActivesRhs<gsMatrix<T> >(rhs));
        return result;
    }

    template <typename T>// to do: rename to getRecipe
    static gsRecipe<T> makePoissonRecipe (const gsFunction<T>  &rhsD, gsSparseMatrix<T> &stiff, gsMatrix<T> &rhs, const gsDofMapper &dofmapper, gsSparseMatrix<T> &rhs_mod, int patch_id)
    {
        gsRecipe<T> result(2);
        result[STIFFNESS].setOperator(new gsGradGradOp<T>());
        result[STIFFNESS].setTestSpace(0);
        result[STIFFNESS].setUnknownSpace(0);
        result[STIFFNESS].setRule(new gsL2GMapped<gsSparseMatrix<T>, gsSparseMatrix<T> >(stiff, rhs_mod, dofmapper,dofmapper,patch_id));

        result[RHS].setOperator(new gsL2TestOp<T>(rhsD));
        result[RHS].setTestSpace(0);
        result[RHS].setUnknownSpace(0);
        result[RHS].setRule(new gsL2GMappedRhs<gsMatrix<T> >(rhs,dofmapper,patch_id));
        return result;
    }

    /**
        \brief provides a Lagrange multiplier for pure Neumann problems
        remember to increase the size of the matrix by one
        here we are assuming a 1 dimensional problem
    **/
    template <typename T>// to do: rename to getZeroAverageMultiplier
    static gsRecipe<T> makeZeroAverageMultiplier (gsSparseMatrix<T> &stiff, const gsDofMapper &dofmapper, int patch_id)
    {
        typedef gsMultiplierWriter<gsSparseMatrix<T>, gsSparseMatrix<T> > MW;
        typedef gsShiftWriter<MW> SW;
        typedef gsDOFMappedMapping<T> DOF;
        typedef gsL2GBase<DOF,gsIdentityMapping, SW> Writer;
        static const gsConstantFunction<T> my_one_func(T(1));
        gsRecipe<T> result(1);
        result[0].setOperator(new gsL2TestOp<T>(my_one_func));
        result[0].setTestSpace(0);
        result[0].setUnknownSpace(0);
        result[0].setRule(new Writer(DOF(dofmapper,patch_id),gsIdentityMapping(), SW(MW(stiff,stiff),0,dofmapper.freeSize())));
        return result;
    }


};

struct gsPoissonBoundaryCookBook { // why not all Poisson related things in one cookbook ?
    // ingredients for Neumann sides
    enum {
        NEUMANN = 0
    };
    // ingredients to use Nitsche method
    enum {
        NITSCHE_NDV  = 0, // normal deriv against value
        NITSCHE_NDD  = 1, // normal deriv against Dirichlet data
        NITSCHE_MASS = 2, // mass matrix on the boundary
        NITSCHE_L2T  = 3  // l2test on the boundary
    };
    // ingredients to compute Dirichlet fixed dofs
    enum  {
        BOUNDARY_MASS = 0,
        BOUNDARY_RHS  = 1
    };
    // ingredients to solve Robin type conditions
    enum {
        ROBIN_MASS = 0,
        ROBIN_RHS  = 1
    };


    /// TODO make the geometry evaluator know about the element size
    /// or the basis evaluator about the support size, so we can
    /// deduce Nitsche constant locally in the operators
    /// this is part of a more general discussion
    template <typename T>
    static gsRecipe<T> makeNitscheRecipe (
            T                     nitsche,
            T                     hparam,
            gsFunction<T>        &dirF,
            const gsDofMapper &dofMap,
            gsSparseMatrix<T>    &sys_m,
            gsMatrix<T>          &sys_rhs,
            unsigned              patch_id,
            boxSide        side_id )
    {
        typedef gsCoeffWriter<gsSparseMatrix<T> >        WM;
        typedef gsCoeffWriter<gsMatrix<T> >              WRHS;
        typedef gsMultiplierWriter<gsSparseMatrix<T> >   WMUL;
        typedef gsCoeffWriter<WMUL>                      MNW;
        typedef gsDOFMappedMapping<T>                    trans;

        gsRecipe<T> result(4);
        result[NITSCHE_MASS].setOperator(new gsBoundaryL2ScalarOp<T>(side_id));
        result[NITSCHE_MASS].setTestSpace(0);
        result[NITSCHE_MASS].setUnknownSpace(0);
        result[NITSCHE_MASS].setRule(new gsL2GBase<trans,trans,WM>(trans(dofMap,patch_id),trans(dofMap,patch_id), WM(sys_m,nitsche/hparam)));

        result[NITSCHE_L2T].setOperator(new gsBoundaryL2TestOp<T>(dirF,side_id));
        result[NITSCHE_L2T].setTestSpace(0);
        result[NITSCHE_L2T].setUnknownSpace(0);
        result[NITSCHE_L2T].setRule(new gsL2GBase<trans,gsIdentityMapping,WRHS>(trans(dofMap,patch_id),gsIdentityMapping(),WRHS(sys_rhs,nitsche/hparam)));

        result[NITSCHE_NDD].setOperator(new gsBoundaryNormalDerTestOp<T>(dirF,side_id));
        result[NITSCHE_NDD].setTestSpace(0);
        result[NITSCHE_NDD].setUnknownSpace(0);
        result[NITSCHE_NDD].setRule(new gsL2GBase<trans,gsIdentityMapping,WRHS>(trans(dofMap,patch_id),gsIdentityMapping(),WRHS(sys_rhs,-1)));

        result[NITSCHE_NDV].setOperator(new gsBoundaryNormalDerValueOp<T>(side_id));
        result[NITSCHE_NDV].setTestSpace(0);
        result[NITSCHE_NDV].setUnknownSpace(0);
        result[NITSCHE_NDV].setRule(new gsL2GBase<trans,trans,MNW>(trans(dofMap,patch_id),trans(dofMap,patch_id), MNW( WMUL(sys_m),-1) ));
        return result;
    }


    template <typename T>
    static gsRecipe<T> makeRobinRecipe (
            gsFunction<T>        &dirF,
            const gsDofMapper &dofMap,
            gsSparseMatrix<T>    &sys_m,
            gsMatrix<T>          &sys_rhs,
            unsigned              patch_id,
            boxSide        side_id )
    {
        typedef gsBaseWriter<gsSparseMatrix<T> >        WM;
        typedef gsBaseWriter<gsMatrix<T> >              WRHS;
        typedef gsDOFMappedMapping<T>                    trans;

        gsRecipe<T> result(2);
        result[ROBIN_MASS].setOperator(new gsBoundaryL2ScalarOp<T>(side_id));
        result[ROBIN_MASS].setTestSpace(0);
        result[ROBIN_MASS].setUnknownSpace(0);
        result[ROBIN_MASS].setRule(new gsL2GBase<trans,trans,WM>(trans(dofMap,patch_id),trans(dofMap,patch_id), WM(sys_m)));

        result[ROBIN_RHS].setOperator(new gsBoundaryL2TestOp<T>(dirF,side_id));
        result[ROBIN_RHS].setTestSpace(0);
        result[ROBIN_RHS].setUnknownSpace(0);
        result[ROBIN_RHS].setRule(new gsL2GBase<trans,gsIdentityMapping,WRHS>(trans(dofMap,patch_id),gsIdentityMapping(),WRHS(sys_rhs)));

        return result;
    }


    template <typename T>
    static gsRecipe<T> makeRecipeForDiricletValues (
            gsFunction<T>        &dirF,
            const gsDofMapper &dofMap,
            gsSparseMatrix<T>    &dirichlet_sys_m,
            gsMatrix<T>          &dirichlet_rhs_m,
            unsigned              patch_id)
    {
        typedef gsBoundaryWriter<gsSparseMatrix<T> > SW;
        typedef gsBoundaryWriter<gsMatrix<T> > SWR;

        typedef gsDOFMappedMapping<T>             DOF;
        typedef gsL2GBase<DOF,DOF,SW>                   MassMapper;
        typedef gsL2GBase<DOF,gsIdentityMapping,SWR>     RhsMapper;

        gsRecipe<T> result(2);
        result[BOUNDARY_MASS].setOperator(new gsBoundaryL2ScalarOp<T>());
        result[BOUNDARY_MASS].setTestSpace(0);
        result[BOUNDARY_MASS].setUnknownSpace(0);
        result[BOUNDARY_MASS].setRule(new MassMapper(DOF(dofMap,patch_id),DOF(dofMap,patch_id),SW(dirichlet_sys_m,dofMap.freeSize(),dofMap.freeSize())));
        result[BOUNDARY_RHS].setOperator(new gsBoundaryL2TestOp<T>(dirF));
        result[BOUNDARY_RHS].setTestSpace(0);
        result[BOUNDARY_RHS].setUnknownSpace(0);
        result[BOUNDARY_RHS].setRule(new RhsMapper(DOF(dofMap,patch_id),gsIdentityMapping() ,SWR(dirichlet_rhs_m,dofMap.freeSize(),0)));
        return result;
    }


    template <typename T>
    static gsRecipe<T> makeNeumannRecipe (
            gsFunction<T>         &rhsF,
            const gsDofMapper  &dofMap,
            gsMatrix<T>           &rhs_m,
            unsigned               patch_id,
            boxSide         side_id )
    {
        // typedef gsShiftWriter<gsSparseMatrix<T> > SW;
        // typedef gsShiftWriter<gsMatrix<T> > SWR;
        
        //typedef gsDOFMappedMapping<T>             DOF;
        //typedef gsL2GBase<DOF,DOF,SW>                   MassMapper;
        //typedef gsL2GBase<DOF,gsIdentityMapping,SWR>     RhsMapper;

        gsRecipe<T> result(1);
        result[NEUMANN].setOperator(new gsBoundaryL2TestOp<T>(rhsF));
        result[NEUMANN].setTestSpace(0);
        result[NEUMANN].setUnknownSpace(0);
        result[NEUMANN].setRule(new gsL2GMappedRhs<gsMatrix<T> >(rhs_m,dofMap,patch_id));
        return result;
    }

};


struct gsPoissonInterfacesCookBook // why not all Poisson related things in one cookbook ?
{
    /// fluxes interfaces for DG patch matching
    /// this is a beer and pizza quest:
    ///  -- change   interface in boxTopology to contain the mapping of the points
    ///  -- provide  iterators for the interface elements: they must iterate in the coarsest grid that is
    ///     a subgrid of both grids
    ///  -- either provide a otherPatch basis evaluator or an ad hoc reimplementation of gsRecipeAssembler
    ///     (the first option being the best)
    /// A maybe working implementation for discontinuous Galerkin coupling is available in PoissonSolver::applyDG

};

} // namespace gismo

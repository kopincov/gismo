/** @file gsTensorAssemblerTest.cpp

    @brief This program tests that the gsRecipeAssembler assembles in 2D really K(x)M+M(x)K
    where M and K are the corresponding 1D matrices.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <iostream>
#include <iomanip>
#include <ctime>


#include <gismo.h>
#include <gismo_dev.h>

#include <gsSolver/gsConjugateGradient.h>
#include <gsSolver/gsKronecker.h>

//#include <gsIO/gsMatrixIO.h>
#include <gsTensor/gsTensorTools.h>

#include <gsAssembler/gsParameterDomainAssembler.h>

using namespace gismo;


int main(int argc, char *argv[])
{

    index_t numRefine = 5;
    index_t degree = 6;
    std::string bcs("nnnn");
    
    gsInfo.precision(3);
    gsInfo << std::fixed;
    
    gsCmdLine cmd("Compares two different poisson assemblers");
    cmd.addInt("r","refine","refinement level", numRefine);
    cmd.addInt("d","degree","polynomial degree", degree);
    cmd.addString("b","bc","boundary conditions", bcs);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Polynomial degree   : " << degree
        << "\nRefinement level    : " << numRefine << std::endl << std::endl;
    
    gsStopwatch time; 
    double t1, t2, t3 = 0.;
    gsGeometry<>::uPtr geo( gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(degree)) );
    gsBasis<>::uPtr basis = geo->basis().clone();
    
    for (int i = 0; i < numRefine; ++i)
        basis->uniformRefine();

    gsFunctionExpr<> f;

    f = gsFunctionExpr<>("(2*pi^2 ) * sin(pi*x) * sin(pi*y)",2);
    
    
    gsMultiPatch<> mp(*geo);
    gsMultiBasis<> mb(*basis);

    gsConstantFunction<> zero(0.0, geo->geoDim());
    
    if ( bcs.length() != 4 )
    {
        gsInfo << "The given boundary conditions are not valid.\n"; return -1;
    }
    for ( int i=0; i<4; ++i )
        if ( bcs[i] != 'd' && bcs[i] != 'n' )
        {
            gsInfo << "The given boundary conditions are not valid.\n"; return -1;
        }
        
    
    gsBoundaryConditions<> bc;
    bc.addCondition( boundary::west,  bcs[0]=='n' ? condition_type::neumann : condition_type::dirichlet, &zero );
    bc.addCondition( boundary::east,  bcs[1]=='n' ? condition_type::neumann : condition_type::dirichlet, &zero );
    bc.addCondition( boundary::south, bcs[2]=='n' ? condition_type::neumann : condition_type::dirichlet, &zero );
    bc.addCondition( boundary::north, bcs[3]=='n' ? condition_type::neumann : condition_type::dirichlet, &zero );
        
    gsInfo << "Assemble in 2D...                              " << std::flush;
    gsOptionList assemblerOptions = gsAssembler<>::defaultOptions();
    gsPoissonAssembler<real_t> assm(mp, mb, bc, f,
                                    (dirichlet::strategy)assemblerOptions.getInt("DirichletStrategy"),
                                    (iFace::strategy)assemblerOptions.getInt("InterfaceStrategy"));
    //assm.computeDirichletDofs();
    assm.assemble();
    gsSparseMatrix<real_t> Kdirect = assm.matrix();
    gsMatrix<real_t> fdirect = assm.rhs();
    
    t1 = time.stop(); gsInfo << "done in " << t1 << " secs.\n";
    
    gsInfo << "Assemble Kroneckerized...                      " << std::flush;

    gsSparseMatrix<> Kkronecker;
    assembleParameterStiffnessForTensorProductSpace(*basis, bc, Kkronecker);
        
    t2 = time.stop(); gsInfo << "done in " << (t2-t1) << " secs.\n";
    
    gsInfo << "Assemble RHS Kroneckerized...                  " << std::flush;

    gsMatrix<> fkronecker;
    assembleParameterMomentsForTensorProduct(*basis, bc, f, fkronecker);
        
    t3 = time.stop(); gsInfo << "done in " << (t3-t2) << " secs.\n";
    
    
    gsSparseMatrix<> Kdiff = Kdirect - Kkronecker;
    gsMatrix<> fdiff = fdirect - fkronecker;
    
    gsInfo << std::endl 
        << "Norm of Kdirect    : " << Kdirect.norm() << std::endl
        << "Norm of Kkronecker : " << Kkronecker.norm() << std::endl
        << "Norm of error      : " << Kdiff.norm() << std::endl
        << "DOFs               : " << Kdirect.rows() << std::endl << std::endl;

    gsInfo << std::endl 
        << "Norm of fdirect    : " << fdirect.norm() << std::endl
        << "Norm of fkronecker : " << fkronecker.norm() << std::endl
        << "Norm of error      : " << fdiff.norm() << std::endl
        << "DOFs               : " << fdirect.rows() << std::endl << std::endl;
    
    return (Kdiff.norm() < 1.e-12 && fdiff.norm() < 1.e-12) ? 0 : -1;
    
}



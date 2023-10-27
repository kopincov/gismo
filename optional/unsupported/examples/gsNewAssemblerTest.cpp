/** @file gsNewAssemblerTest.cpp

    @brief Compares several assembler:xs

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

#include <gsAssembler/gsParameterDomainAssembler.h>

#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>


//#include <gsIO/gsMatrixIO.h>
#include <gsTensor/gsTensorTools.h>

using namespace gismo;

class gsRecipeAssemblerGenericSecondOrderOp : public gsRecipeAssemblerPoisson
{
protected:
    gsFunction<real_t> *m_A;
    gsFunction<real_t> *m_b;
    gsFunction<real_t> *m_c;
public:
     gsRecipeAssemblerGenericSecondOrderOp( const gsPoissonPde<real_t> &pde, gsFunction<real_t> *A=NULL,gsFunction<real_t> *b=NULL,gsFunction<real_t> *c=NULL)
        : gsRecipeAssemblerPoisson(pde), m_A(A), m_b(b), m_c(c)
    {}

protected:
    virtual gsRecipe<real_t>    getPatchRecipe     (index_t patch)
    {
        gsRecipe<real_t>            result;
        gsRecipeIngredient<real_t>  ingr;

        ingr.setOperator(new gsGenericSecondOrderOp<real_t>(m_A,m_b,m_c));
        //ingr.setOperator(new gsGradGradOp<real_t>());
        ingr.setTestSpace(0);
        ingr.setUnknownSpace(0);
        ingr.setRule(getSysWriter());
        result.add(ingr);

        ingr.setOperator(new gsL2TestOp<real_t>(*m_pde.rhs()));
        ingr.setTestSpace(0);
        ingr.setUnknownSpace(0);
        ingr.setRule(getRhsWriter());
        result.add(ingr);

        return result;

    }

};


void bcChoose( char bc, condition_type::type& bc_type, gsFunction<> *& bc_func, gsFunction<> * bc_func_dirichlet, gsFunction<> * bc_func_neumann )
{
    if( bc == 'd' )
    {
        bc_type = condition_type::dirichlet;
        bc_func = bc_func_dirichlet;
    }
    else if( bc == 'n' )
    {
        bc_type = condition_type::neumann;
        bc_func = bc_func_neumann;
    }
    else
    {
        gsInfo << "Invalid boundary condition. Allowed are: dirichlet (d), neumann (n) and mixed (dd, dn, nd, nn, ddd, ddn, dnd, dnn, ndd, ndn, nnd, nnn).\n";
        exit(-1);
    }
}

int main(int argc, char *argv[])
{

    srand( (unsigned)time( NULL ) );

    index_t geoIndex = 1;
    std::string boundaryCondition("n");
    index_t numRefine = 2;
    index_t degree = 7;
    std::string show1 = "K1";
    std::string show2 = "K4"; 

    gsCmdLine cmd("Compares 4 different assemblers.");
    cmd.addInt("g", "geometry", "Specification of the geometry (overrides dimension)", geoIndex);
    cmd.addString("b", "boundary-condition", "Boundary condition", boundaryCondition);
    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("p", "degree", "Degree of the B-spline discretization space", degree);
    cmd.addString("", "show1", "First matrix to be shown", show1);
    cmd.addString("", "show2", "Second matrix to be shown", show2);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    /******************** Define Geometry ********************/

    gsStopwatch time;    
    
    gsGeometry<>::uPtr geo;

    switch ( geoIndex )
    {
        case 1: geo = gsNurbsCreator<>::BSplineUnitInterval(static_cast<short_t>(degree));  break;
        case 2: geo = gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(degree));     break;
        case 3: geo = gsNurbsCreator<>::BSplineCube(static_cast<short_t>(degree));          break;
        //case 4: geo = approximateQuarterAnnulus(degree);              break;
        //case 5: geo = BSplineMySquare(degree);                        break;
        default: gsWarn << "Invalid geometry. Allowed are:\n"
            << "1: unit interval\n"
            << "2: unit square\n"
            << "3: unit cube\n";
            //<< "4: approximate quarter annulus\n"
            //<< "5: unit square with one additional refinement in x-direction\n";
        return -1;
    }

    gsFunctionExpr<> f, g;

    switch (geo->geoDim())
    {
        case 1:
            f = gsFunctionExpr<>("(pi^2 ) * sin(pi*x)",1);
            g = gsFunctionExpr<>("sin(pi*x)",1);
            break;
        case 2:
            f = gsFunctionExpr<>("(2*pi^2 ) * sin(pi*x) * sin(pi*y)",2);
            g = gsFunctionExpr<>("sin(pi*x) * sin(pi*y)",2);
            break;
        case 3:
            f = gsFunctionExpr<>("(3*pi^2 ) * sin(pi*x) * sin(pi*y) * sin(pi*z)",3);
            g = gsFunctionExpr<>("sin(pi*x) * sin(pi*y) * sin(pi*z)",3);
            break;
        default:
            gsWarn << "Invalid geometry dimension.\n";
            return -1;
    }

    if (numRefine < 0)
    {
        gsWarn << "Number of refinements must not be negative.\n"; return -1;
    }
    
    // set up boundary conditions

    gsConstantFunction<> zero(0.0, geo->geoDim());
    //gsConstantFunction<> one (1.0, geo->geoDim());

    gsBoundaryConditions<> bc;
    condition_type::type bc_type;
    gsFunction<> * bc_func;
    
    if (boundaryCondition.length() == 1)
        boundaryCondition = std::string(geo->geoDim(), boundaryCondition[0]);

    if ( (index_t)boundaryCondition.length() != geo->geoDim() )
        boundaryCondition = "x"; // Let the bcChoose do the work
    
    bcChoose( boundaryCondition[0], bc_type, bc_func, &g, &zero );
    bc.addCondition( boundary::west,  bc_type, bc_func );
    bc.addCondition( boundary::east,  bc_type, bc_func );
    if (geo->geoDim() >= 2)
    {
        bcChoose( boundaryCondition[1], bc_type, bc_func, &g, &zero );
        bc.addCondition( boundary::south, bc_type, bc_func );
        bc.addCondition( boundary::north, bc_type, bc_func );
    }
    if (geo->geoDim() >= 3)
    {
        bcChoose( boundaryCondition[2], bc_type, bc_func, &g, &zero );
        bc.addCondition( boundary::front, bc_type, bc_func );
        bc.addCondition( boundary::back,  bc_type, bc_func );
    }

    gsBasis<> * tbasis = geo->basis().clone().release();
    
    for (int i = 0; i < numRefine; ++i)
       tbasis->uniformRefine();

    gsMultiPatch<> mp(*geo);
    gsMultiBasis<> mb(*tbasis);
    gsGenericAssembler<real_t> genassm(mp, mb, gsAssembler<>::defaultOptions(), &bc);
    
    gsInfo << "Assembling stiffness matrix K1 with gsGenericAssembler... " << std::flush;
    genassm.assembleStiffness();
    gsSparseMatrix<real_t> K1 = genassm.fullMatrix();
    gsInfo << "done, " << K1.rows() << " dofs." << std::endl;
    
    gsInfo << "Assembling mass matrix M1 with gsGenericAssembler... " << std::flush;
    genassm.assembleMass();
    gsSparseMatrix<real_t> M1 = genassm.fullMatrix();
    gsInfo << "done, " << M1.rows() << " dofs." << std::endl;

    gsInfo << "Assemble stiffness matrix K2 with gsPoissonAssembler..." << std::flush;
    gsPoissonAssembler<> passm(mp, mb, bc, f);
    passm.assemble();
    gsSparseMatrix<real_t> K2 = passm.matrix();
    gsInfo << "done, " << K2.rows() << " dofs." << std::endl;


    gsInfo << "Assembling stiffness matrix K3 with gsRecipeAssemblerPoisson... " << std::flush;
    gsPoissonPde<real_t> pde(mp,bc,f,&g);
    gsRecipeAssemblerPoisson assembler(pde);
    assembler.setDirichletStrategy(dirichlet::elimination);
    assembler.setZeroAverage(false);
    std::vector<gsPhysicalSpace*> phySpace;
    phySpace.push_back(new gsPhysicalSpaceScalar(mb,mp,INVERSE_COMPOSITION));
    assembler.setSpace(phySpace);
    assembler.assemble();
    gsSparseMatrix<real_t> K3 = assembler.getSystemMatrix();
    gsInfo << "done, " << K3.rows() << " dofs." << std::endl;

    /*
    gsInfo << "Assembling mass matrix M3 with gsRecipeAssemblerGenericSecondOrderOp... " << std::flush;
    gsFunction<>* c= new gsConstantFunction<>( 1., 2 );
    gsRecipeAssemblerGenericSecondOrderOp assemblerM(pde, NULL, NULL, c);
    assemblerM.setDirichletStrategy(dirichlet::elimination);
    assemblerM.setZeroAverage(false);
    assemblerM.setSpace(phySpace);
    assemblerM.assemble();
    gsSparseMatrix<real_t> M3 = assemblerM.getSystemMatrix();
    gsInfo << "done, " << M3.rows() << " dofs." << std::endl;
    //*/

    gsInfo << "Assembling stiffness matrix K4 with assembleParameterStiffness..." << std::flush;
    gsSparseMatrix<real_t> K4;
    assembleParameterStiffness( *tbasis, K4 );
    gsInfo << "done, " << K4.rows() << " dofs." << std::endl;


    gsInfo << "Assembling mass matrix M4 with assembleParameterMass..." << std::flush;
    gsSparseMatrix<real_t> M4;
    assembleParameterMass( *tbasis, M4 );
    gsInfo << "done, " << M4.rows() << " dofs." << std::endl;


    real_t x; bool ok = true;

    gsInfo << "Error K1-K2: " << (x=(K1-K2).norm()) << std::endl; ok &= (x<1.e-5); 
    gsInfo << "Error K1-K3: " << (x=(K1-K3).norm()) << std::endl; ok &= (x<1.e-5);
    gsInfo << "Error K1-K4: " << (x=(K1-K4).norm()) << std::endl; ok &= (x<1.e-5);
    gsInfo << "Error K2-K3: " << (x=(K2-K3).norm()) << std::endl; ok &= (x<1.e-5);
    gsInfo << "Error K2-K4: " << (x=(K2-K4).norm()) << std::endl; ok &= (x<1.e-5);
    gsInfo << "Error K3-K4: " << (x=(K3-K4).norm()) << std::endl; ok &= (x<1.e-5);

    gsInfo << std::endl;

    gsInfo << "Error M1-M4: " << (x=(M1-M4).norm()) << std::endl; ok &= (x<1.e-5);
    //gsInfo << "Error M1-M3: " << (x=(M1-M3).norm()) << std::endl; ok &= (x<1.e-5);
    //gsInfo << "Error M3-M4: " << (x=(M3-M4).norm()) << std::endl; ok &= (x<1.e-5);


    gsInfo << std::endl;


    gsSparseMatrix<real_t> res1, res2;

    if ( show1 == "K1" )
       res1=K1;
    else if ( show1 == "K2" )
       res1=K2;
    else if ( show1 == "K3" )
       res1=K3;
    else if ( show1 == "K4" )
       res1=K4;
    else if ( show1 == "M1" )
       res1=M1;
    //else if ( show1 == "M3" )
    //   res1=M3;
    else if ( show1 == "M4" )
       res1=M4;
    else if ( show1 != "" )
    {
       gsInfo << "show1 has an invalid value" << std::endl; return -1;
    }

    if ( show1.length() > 0 )
	gsInfo << show1 << std::endl << res1 << std::endl;

    if ( show2 == "K1" )
       res2=K1;
    else if ( show2 == "K2" )
       res2=K2;
    else if ( show2 == "K3" )
       res2=K3;
    else if ( show2 == "K4" )
       res2=K4;
    else if ( show2 == "M1" )
       res2=M1;
    //else if ( show2 == "M3" )
    //   res2=M3;
    else if ( show2 == "M4" )
       res2=M4;
    else if ( show2 != "" ) 
    {
       gsInfo << "show2 has an invalid value" << std::endl; return -1;
    }

    if ( show2.length() > 0 )
        gsInfo << show2 << std::endl << res2 << std::endl;
  
    if ( show1.length() > 0 && show2.length() > 0 )
        gsInfo << show1 << "-" << show2 << std::endl << ( res1 - res2 ) << std::endl;

    delete tbasis;
    freeAll( phySpace );
    
    if ( ok )
    {
        gsInfo << "Tests passed" << std::endl;
        return 0;
    }
    else
    {
        gsInfo << "Tests failed" << std::endl;
        return -1;
    }

}

 

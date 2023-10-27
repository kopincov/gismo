
/** @file biharmonicModified.cpp

    @brief A modified Biharmonic example.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

# include <gismo.h>
# include <gsAssembler/gsBiharmonicAssembler.h>

using namespace gismo;

int main(int argc, char *argv[])
{

    index_t numRefine = 0;
    if (argc >= 2)
        numRefine = atoi(argv[1]);

    gsFunctionExpr<> sol44("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsFunctionExpr<> f44("cos(4*pi*x)*cos(4*pi*y) + 64*pi*pi*cos(4*pi*x)*cos(4*pi*y) + 1024*pi^4*cos(4*pi*x)*cos(4*pi*y) - 256*pi^4*cos(4*pi*x) - 32*pi*pi*cos(4*pi*x) - cos(4*pi*x) - 256*pi^4*cos(4*pi*y) - 32*pi*pi*cos(4*pi*y) - cos(4*pi*y) + 1",2);
    gsFunctionExpr<> laplace44 ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);

    gsMultiPatch<> geo(*gsNurbsCreator<>::BSplineSquareDeg(static_cast<short_t>(2)));
    //geo( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    //geo.uniformRefine();
    //geo.uniformRefine();

    /*
    gsMultiBasis<> basisTmp = gsMultiBasis<>(geo);
    gsTensorBSplineBasis<2,real_t> * tbasis;
    tbasis= basis.basis(0).clone());

    gsInfo << "test\n" << tbasis->knots(0) << "\n";
    gsInfo << "Coefs:\n"<<geo.patch(0).coefs() << "\n";

    gsMatrix<> newCoefs = geo.patch(0).coefs();
    //outher coefs
    newCoefs(1,0) = 0; newCoefs(4,0) = 1;

    newCoefs(7,0) = 0; newCoefs(10,0) = 1;
    newCoefs(6,1) = 0; newCoefs(7,1) = 0; newCoefs(8,1) = 0; newCoefs(9,1) = 0; newCoefs(10,1) = 0; newCoefs(11,1) = 0;

    newCoefs(13,0) = 0; newCoefs(16,0) = 1;

    newCoefs(19,0) = 0; newCoefs(22,0) = 1;

    newCoefs(25,0) = 0; newCoefs(28,0) = 1;
    newCoefs(24,1) = 1; newCoefs(25,1) = 1; newCoefs(26,1) = 1; newCoefs(27,1) = 1; newCoefs(28,1) = 1; newCoefs(29,1) = 1;

    newCoefs(31,0) = 0; newCoefs(34,0) = 1;
    //Inner  coefs
    real_t a = 0.25;
    real_t b = 0.75;
    newCoefs(2,0) = a; newCoefs(3,0) = b;

    newCoefs(8,0) = a; newCoefs(9,0) = b;

    newCoefs(14,0) = a; newCoefs(15,0) = b;
    newCoefs(12,1) = a; newCoefs(13,1) = a; newCoefs(14,1) = a; newCoefs(15,1) = a; newCoefs(16,1) = a; newCoefs(17,1) = a;

    newCoefs(20,0) = a; newCoefs(21,0) = b;
    newCoefs(18,1) = b; newCoefs(19,1) = b; newCoefs(20,1) = b; newCoefs(21,1) = b; newCoefs(22,1) = b; newCoefs(23,1) = b;

    newCoefs(26,0) = a; newCoefs(27,0) = b;

    newCoefs(32,0) = a; newCoefs(33,0) = b;

    gsInfo << "Coefs:\n"<<newCoefs << "\n";
    gsKnotVector<real_t> KV = tbasis->knots(0);
    gsTensorBSpline<2,real_t> * newSpline = new gsTensorBSpline<2,real_t>(KV,KV, give(newCoefs));
    gsMultiPatch<> newGeo( *newSpline );
    */

    //gsMultiBasis<> basis = gsMultiBasis<>(newGeo);
    gsMultiBasis<> basis = gsMultiBasis<>(geo);
    //gsMultiPatch<> geoAn( *gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    //gsMultiBasis<> basis = gsMultiBasis<>(geoAn);

    for (int i = 0; i < numRefine; ++i)
    {
        basis.uniformRefine();
    }

    gsBoundaryConditions<> bcInfo;
    gsBoundaryConditions<> bcInfo2;
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &sol44);//Annulus: small arch lenght
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &sol44);//Annulus: Large arch lenght
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &sol44);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &sol44);
    bcInfo2.addCondition( boundary::west,  condition_type::neumann, &laplace44);//Annulus: small arch lenght
    bcInfo2.addCondition( boundary::east,  condition_type::neumann, &laplace44);//Annulus: Large arch lenght
    bcInfo2.addCondition( boundary::north, condition_type::neumann, &laplace44);
    bcInfo2.addCondition( boundary::south, condition_type::neumann, &laplace44);


    /////////////////// Setup solver ///////////////////
    //Initilize Solver
    gsBiharmonicAssembler<real_t> BiharmonicAssembler(geo,basis,bcInfo,bcInfo2,f44,
                                                       dirichlet::elimination, iFace::dg);


    gsInfo<<"Assembling..." << "\n";
    BiharmonicAssembler.assemble();

    gsInfo << "Number of DOFs: " << BiharmonicAssembler.numDofs() << "\n";


    // Initialize the congugate gradient solver
    gsInfo<<"Solving...\n";

    /*
    gsSparseSolver<>::QR solver;
    solver.analyzePattern(BiharmonicAssembler.matrix() );
    solver.factorize(BiharmonicAssembler.matrix());
    gsDebug << "   Eigen  lastErrorMessage: " << solver.lastErrorMessage () << "\n";
    gsMatrix<> solVector= solver.solve(BiharmonicAssembler.rhs());

    gsSparseSolver<>::LU solverLU;
    solverLU.analyzePattern(BiharmonicAssembler.matrix());
    solverLU.factorize(BiharmonicAssembler.matrix());
    //Use the factors to solve the linear system
    gsMatrix<> solVector = solverLU.solve(BiharmonicAssembler.rhs());
    */


    gsConjugateGradient<> solver(BiharmonicAssembler.matrix()); // without precond.
    solver.setMaxIterations(20000);
    solver.setTolerance(1e-12);
    gsMatrix<> solVector;
    solVector.setZero(BiharmonicAssembler.numDofs(),1);
    solver.solve(BiharmonicAssembler.rhs(),solVector);

    gsInfo <<"Number of interations :" <<solver.iterations() << "\n";

    gsMultiPatch<> mpsol;
    BiharmonicAssembler.constructSolution(solVector, mpsol);
    gsField<> solF(BiharmonicAssembler.patches(), mpsol);

    real_t l2error = solF.distanceL2(sol44);
    gsInfo << "The L2 error of the solution: " << l2error << "\n";

    // Plot solution in paraview
    bool plot = true;
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>(solF, "Biharmonic2d", 1000);
        const gsField<> exact( geo, sol44, false );
        gsWriteParaview<>( exact, "Biharmonic2d_exact", 1000);

    }

    gsInfo << "Test is done: Cleaning up..." << "\n"; //freeAll(m_bconditions);

    gsInfo << "Test is done: Exiting" << "\n";

    return  0;
}


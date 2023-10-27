#include <iostream>

#include <gismo.h>


#include <gsAssembler/gsStokesAssembler.h>
#include <gsAssembler/gsPdeAssembler.h>
#include <gsSolver/gsSolverUtils.h>


;
using namespace gismo;

int main(int argc, char *argv[])
{

    int numRefine = 1;
    if (argc == 2)
        numRefine = atoi(argv[1]);

    bool plot = false; // If set to true, paraview file is generated and launched on exit
    /*try
    {
    gsCmdLine cmd("Tutorial 01 shows the use of BSpline curves.");
    gsArgSwitch ap("", "plot", "Plot result in ParaView format", cmd);
    cmd.parse(argc,argv);
    plot       = ap.getValue();
    } catch ( gsArgException& e )
    { gsInfo << "Error: " << e.error() << " " << e.argId() << "\n"; return -1; }
    */

    // Source function
    //gsFunctionExpr<> f("2*pi^2*sin(pi*x)*sin(pi*y)",3) ;
    //gsFunctionExpr<> f(" 2*pi^3*sin(2*pi*y)*(2*(sin(pi*x)^2)- cos(2*pi*x)) + 2*pi*cos(2*pi*x)",
    //                          "-2*pi^3*sin(2*pi*x)*(2*(sin(pi*y)^2)- cos(2*pi*y))",3) ;
    //gsFunctionExpr<> f("2*cos(0.0)", "0.0", "0.0",3) ;
    gsFunctionExpr<> f("2*4*pi*pi*sin(2*pi*x)*cos(2*pi*y) + 4*pi*pi*sin(2*pi*x)",
                        "-2*4*pi*pi*cos(2*pi*x)*sin(2*pi*y)", "0.0",3) ;

    // Exact solution
    //gsFunctionExpr<> g("pi*sin(2*pi*y)*(sin(pi*x)^2)", "-pi*sin(2*pi*x)*(sin(pi*y)^2)",2);
    //gsFunctionExpr<> g("y*(1-y)", "0", "0",3);
    gsFunctionExpr<> g("sin(2*pi*x)*cos(2*pi*y)", "-cos(2*pi*x)*sin(2*pi*y)", "0.0",3);


    //gsFunctionExpr<> p0("2*pi*cos(2*pi*x)",3);

    gsInfo<<"Source function "<< f << "\n";
    gsInfo<<"Exact solution "<< g <<".\n" << "\n";


    // Define Geometry
    gsMultiPatch<> * patches;
    //patches = new gsMultiPatch<>(*gsNurbsCreator<>::BSplineFatQuarterAnnulus() );
    //patches = new gsMultiPatch<>(*gsNurbsCreator<>::BSplineSquareDeg(2));
    // Unit cube
    patches = new gsMultiPatch<>(*gsNurbsCreator<>::BSplineCube(static_cast<short_t>(1)));
    // Define Geometry (Unit square with 4 patches)
    //patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);

    std::vector< gsMultiBasis<>* >  refine_bases;
    refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for velocity_x
    refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for velocity_y
    refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for velocity_z
    refine_bases.push_back(new gsMultiBasis<>( *patches ));//Basis for pressure

    // Define discretization space by refining the basis of the geometry
    for (int i = 0; i < numRefine; ++i)
    {
        refine_bases[0]->uniformRefine();
        refine_bases[1]->uniformRefine();
        refine_bases[2]->uniformRefine();
        refine_bases[3]->uniformRefine();
    }
    gsDebug <<  " refine_bases.size() " << refine_bases.size() << "\n";

    refine_bases[0]->degreeElevate();
    refine_bases[1]->degreeElevate();
    refine_bases[2]->degreeElevate();
    //refine_bases[0][i]->degreeElevate();

    //refine_bases[1][i]->degreeElevate();
    //refine_bases[1][i]->degreeReduce();



    gsInfo << "Degree of basis for velocity: "<< refine_bases[0]->degree()<< "\n";
    gsInfo << "Degree of basis for velocity: "<< refine_bases[1]->degree()<< "\n";
    gsInfo << "Degree of basis for velocity: "<< refine_bases[2]->degree()<< "\n";
    gsInfo << "Degree of basis for pressure: "<< refine_bases[3]->degree()<< "\n";
    gsInfo << "basis u" << *refine_bases[0] << "\n";
    gsInfo << "basis p" << *refine_bases[1] << "\n";


    gsInfo << "Number of patches is " << patches->nPatches() << "\n";


    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;
    //gsBoundaryConditions<> bcInfo_p;

    //Dirichlet BCs

    gsFunctionExpr<> h_neu("-2*pi*(1+cos(2*pi*y))","0.0", "0.0",3);

    //front back
    bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &g);
    bcInfo.addCondition( boundary::east,  condition_type::dirichlet, &g);
    bcInfo.addCondition( boundary::north, condition_type::dirichlet, &g);
    bcInfo.addCondition( boundary::south, condition_type::dirichlet, &g);
    bcInfo.addCondition( boundary::front, condition_type::dirichlet, &g);
    bcInfo.addCondition( boundary::back, condition_type::dirichlet, &g);

    //bcInfo.addCondition( boundary::west,  condition_type::neumann, &h_neu);
    //bcInfo.addCondition( boundary::west,  condition_type::dirichlet, &g);


    //std::vector< gsMultiBasis<> >refine_Mbases ;

    std::vector<gsBoundaryConditions<>* > vecBcInfo;
    vecBcInfo.resize(2);
    vecBcInfo[0] = &bcInfo;
    gsBoundaryConditions<>*  prt_null = NULL;
    vecBcInfo[1] = prt_null;//&bcInfo_p;


    /////////////////// Setup solver ///////////////////
    //Initilize Solver

    // Test contructor 1, 2, 3 and 4
    gsDebug << "-------------Testing contructors--------------\n";
    //gsStokesAssembler<> StokesSolver1(*patches, bcInfo, f);
    //gsDebug << "Contructor 1: OK!\n" ;
    //gsStokesAssembler<> StokesSolver2(*patches, vecBcInfo, f);
    //gsDebug << "Contructor 2: OK!\n" ;
    //gsStokesAssembler<> StokesSolver3(*patches, bcInfo,refine_bases, f);
    //gsDebug << "Contructor 3: OK!\n" ;
    gsStokesAssembler<> StokesSolver4(*patches, vecBcInfo, refine_bases, f);
    gsDebug << "Contructor 4: OK!\n" ;

    gsDebug << "-------------Testing << operator--------------\n";
    //gsDebug << "\n" << "Print Stokes solver 1:\n" <<StokesSolver1<< "\n";
    //gsDebug << "\n" << "Print Stokes solver 2:\n" << StokesSolver2<< "\n";
    //gsDebug << "\n" << "Print Stokes solver 3:\n" << StokesSolver3<< "\n";
    gsDebug << "\n" << "Print Stokes solver 4:\n" << StokesSolver4<< "\n";


    StokesSolver4.setDirichletStrategy(dirichlet::elimination);
    //StokesSolver4.setDirichletStrategy(dirichlet::nitsche);
    StokesSolver4.setInterfaceStrategy(iFace::glue);

    gsDebug << "-------------Testing initialize()--------------\n";
    //StokesSolver1.initialize();
    //gsDebug << "\n\n---------------------------\n\n";
    //StokesSolver2.initialize();
    //gsDebug << "\n\n---------------------------\n\n";
    //StokesSolver3.initialize();
    //gsDebug << "\n\n---------------------------\n\n";
    StokesSolver4.initialize();
    gsDebug << "Testing initializing done!\n";

    gsDebug << "-------------Testing assemble()--------------\n";
    StokesSolver4.assemble();

    gsDebug << "-------------Testing solveSystem()--------------\n";
    StokesSolver4.solveSystem();

    gsDebug << "-------------Testing recontructionSystem()--------------\n";
    StokesSolver4.reconstructSolution();


    gsDebug << "-------------Testing recontructionSystem()--------------\n";
    // Access the solutions
    const gsField<> & sol_u = StokesSolver4.solution(0);
    const gsField<> & sol_p = StokesSolver4.solution(1);


    //double h = math::pow( (double) refine_bases[0]->size(), -1.0 / refine_bases[0]->dim() );

    // Plot solution in paraview
    int result = 0;
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( sol_u, "Stokes2dvelocity", 10000);
        gsWriteParaview<>( sol_p, "Stokes2dpressure", 10000);
        //gsField<> exact( StokesSolver4.patches(), g, false );
        gsField<> exact( *patches , g , false ) ;
        gsWriteParaview<>( exact, "Stokes2d_exact", 10000);

        // Run paraview
        result = system("paraview Stokes2dpressure.pvd &");
    }

    //double conditionnumber =  gsSolverUtils<>::conditionNumber(StokesSolver4.systemMatrix());
    //gsInfo << "the condotion number is : " << conditionnumber << "\n";

    freeAll( refine_bases );
    delete patches;

    gsInfo << "Test New sol" << "\n";
    gsInfo << "Test is done: Exiting" << "\n";
    return  result;

}

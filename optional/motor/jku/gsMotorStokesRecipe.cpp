#include <iostream>

#include <gismo.h>
#include <gismo_dev.h>

#include <gsCore/gsField.h>
#include <gsCore/gsFunctionExpr.h>

#include <gsRecipeAssembler/gsRecipeAssemblerStokes.h>
#include <gsRecipeAssembler/gsRecipeAssemblerDistance.h>

#include "gsMultipatchFitting.h"

using namespace gismo;

gsMatrix<> solve(const gsSparseMatrix<> &sys, const gsMatrix<> &rhs)
{
    gsSparseSolver<>::QR  solver;
    solver.analyzePattern( sys );
    solver.factorize     ( sys );
    return solver.solve( rhs );
}

gsBoxTopology getTopology(const gsMultiPatch<>& mp)
{

    std::vector< patchSide > boundaries = mp.boundaries();
    std::vector< boundaryInterface > interfaces = mp.interfaces();
    
    gsBoxTopology topology(mp.dim(), mp.size(), boundaries, interfaces);

    return topology;
}

std::vector< gsBasis<>* > copyComponent(std::vector< std::vector<gsBasis<>*> >& bases,
                                        const int component)
{
    std::vector< gsBasis<>* > result;
    for (std::size_t i = 0; i != bases[component].size(); i++)
    {
        result.push_back(bases[component][i]->clone().release());
    }

    return result;
}


gsMultiBasis<> getVelocityBasis(std::vector< std::vector<gsBasis<>*> >& bases,
                                           gsBoxTopology& topology)
{
    std::vector< gsBasis<>* > veloc_basis_vector = copyComponent(bases, 0);
    gsMultiBasis<> veloc_basis(veloc_basis_vector, topology);
    //repairInterfaces(veloc_basis); // this doesn't do anything, it is here just to be sure

    return veloc_basis;
}

gsMultiBasis<> getPressureBasis(std::vector< std::vector<gsBasis<>*> >& bases,
                                gsBoxTopology& topology)
{
    std::vector< gsBasis<>* > pressure_basis_vector = copyComponent(bases, 1);
    gsMultiBasis<> press_basis(pressure_basis_vector, topology);
    //repairInterfaces(press_basis);  // this doesn't do anything, it is here just to be sure

    return press_basis;
}

// ================================================================================
// reconstruction version 2

gsMultiPatch<> buildSolution_2(gsMultiBasis<>& basis,
                             const gsBoundaryConditions<>& bc,
                             const gsMatrix<>& coefs)
{
    gsDofMapper dofMapper;
    // basis.getMapper(true, bc, 0, dofMapper);
    basis.getMapper(true, dofMapper);

    dofMapper.print();
    
    std::vector< gsGeometry<>* > patches;
    
    for (std::size_t patch = 0; patch != basis.nBases(); patch++)
    {
        gsBasis<>::uPtr patchBasis = basis[patch].clone();

        const index_t size = patchBasis->size();
        gsMatrix<> patchCoefs(size, coefs.cols());
        patchCoefs.setConstant(100);
        
        for (index_t i = 0; i != size; i++)
        {
            const index_t globalI = dofMapper.index(i, static_cast<index_t>(patch));
            patchCoefs.row(i) = coefs.row(globalI);
        }

        gsGeometry<>::uPtr geom = patchBasis->makeGeometry(patchCoefs);
        patches.push_back(geom.release());
    }

    return gsMultiPatch<>(patches);
}

gsMultiPatch<> buildVelocitySolution_2(gsMultiBasis<>& basis,
                                       const gsBoundaryConditions<>& bc,
                                       const gsMatrix<>& coefs)
{

    const index_t size = coefs.rows() / 2;
    gsMatrix<> newCoefs(size, 2);
    newCoefs.block(0, 0, size, 1) = coefs.block(0, 0, size, 1);
    newCoefs.block(0, 1, size, 1) = coefs.block(size, 0, size, 1);

    return buildSolution_2(basis, bc, newCoefs);
}

gsMultiPatch<> buildPressureSolution_2(gsMultiBasis<>& basis,
                                       const gsBoundaryConditions<>& bc,
                                       const gsMatrix<>& coefs)
{
    return buildSolution_2(basis, bc, coefs);
}

// ================================================================================

// ================================================================================
// reconstruction version 1

gsMultiPatch<> buildSolution(gsMultiBasis<>& basis,
                             const gsMatrix<>& coefs)
{
    std::vector< gsGeometry<>* > patches;
    int start = 0;
    
    for (std::size_t patch = 0; patch != basis.nBases(); patch++)
    {
        gsBasis<>::uPtr patchBasis = basis[patch].clone();
        
        const index_t size = patchBasis->size();
        gsMatrix<> patchCoefs(size, coefs.cols());
        patchCoefs.setZero();

        patchCoefs.block(0, 0, size, coefs.cols()) =
            coefs.block(start, 0, size, coefs.cols());

        start += size;

        gsGeometry<>::uPtr geom = patchBasis->makeGeometry(patchCoefs);
        patches.push_back(geom.release());
    }

    return gsMultiPatch<>(patches);
}
                             
gsMultiPatch<> buildVelocitySolution(gsMultiBasis<>& basis,
                                     const gsMatrix<>& coefs,
                                     std::vector<gsPhysicalSpace*>& phySpace)
{
    gsMapper* veloc_mapper = phySpace[0]->getMapper();

    gsMatrix<> veloc_result;
    veloc_mapper->mapToSourceCoefs(coefs, veloc_result);

    const index_t size = veloc_result.rows() / 2;
    gsMatrix<> newCoefs(size, 2);
    newCoefs.block(0, 0, size, 1) = veloc_result.block(0, 0, size, 1);
    newCoefs.block(0, 1, size, 1) = veloc_result.block(size, 0, size, 1);
    
    return buildSolution(basis, newCoefs);
}

gsMultiPatch<> buildPressureSolution(gsMultiBasis<>& basis,
                                     const gsMatrix<>& coefs,
                                     std::vector<gsPhysicalSpace*>& phySpace)
{
    gsMapper* press_mapper = phySpace[1]->getMapper();

    gsMatrix<> press_result;
    press_mapper->mapToSourceCoefs(coefs, press_result);
    
    return buildSolution(basis, press_result);
}

// ================================================================================


void uniformRefine(gsBasis<>* basis)
{
    basis->uniformRefine();
}


int main(int argc, char *argv[])
{
    std::string multipatchFile(MOTOR_DATA_DIR "jku/airPassageParameterization.xml");
    std::string output("stokes");
    
    // Number for h-refinement of the computational (trial/test) basis.
    int numRefine  = 1;

    // Number for p-refinement of the computational (trial/test) basis.
    int numElevate = 1;

    gsCmdLine cmd("Solve the Stokes equation");
    cmd.addString("m", "m", "Multipatch file", multipatchFile);
    cmd.addInt("", "numElevate", "Number for p-refinement of the computational (trial/test) basis", numElevate);
    cmd.addInt("", "numRefine", "Number for h-refinement of the computational (trial/test) basis", numRefine);
    cmd.addString("o", "out", "Output file", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "---------------------------------------------------------\n\n"
           << "Input Arguments: \n\n"
           << "Multipatch file: " << multipatchFile << "\n"
           << "numRefine:       " << numRefine << "\n"
           << "numElevate:      " << numElevate << "\n"
           << "Output:          " << output << "\n"
           << "---------------------------------------------------------\n"
           << std::endl;

    // Source function
    gsFunctionExpr<> f("0.25*2*4*pi*pi*sin(2*pi*x)*cos(2*pi*y) + 4*pi*pi*sin(2*pi*x)",
                              "-0.25*2*4*pi*pi*cos(2*pi*x)*sin(2*pi*y)", 2) ;

    // Exact solution
    gsFunctionExpr<> g("sin(2*pi*x)*cos(2*pi*y)+pi/10", "-cos(2*pi*x)*sin(2*pi*y)-pi/10", 2);
    gsFunctionExpr<> p0("2*pi*cos(2*pi*x)", 1);

    gsInfo<<"Source function: "<< f << "\n";
    gsInfo<<"Exact solution: " << g << "\n"<< "\n";

/*
    // Geometry case: Out_7patches
    // Define Geometry
    gsFileData<> fd(multipatchFile);
    gsMultiPatch<>* patches = fd.getAnyFirst< gsMultiPatch<> >().release();
    patches->computeTopology();

    gsBoundaryConditions<> bcInfo;

    bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(1, boundary::east,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(2, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(3, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(3, boundary::east,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(4, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(4, boundary::north, condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(5, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(5, boundary::east,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(6, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(6, boundary::south, condition_type::dirichlet, &g, 0);
*/

    // Geometry case: Out_4patches
    // Define Geometry
    gsFileData<> fd(multipatchFile);
    gsMultiPatch<>* patches = fd.getAnyFirst< gsMultiPatch<> >().release();
    patches->computeTopology();
    gsInfo << "Details of geometry:\n" << patches->detail() << std::endl;

    gsBoundaryConditions<> bcInfo;

    bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(1, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(1, boundary::east,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(2, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(2, boundary::east,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(3, boundary::west,  condition_type::dirichlet, &g, 0);
    bcInfo.addCondition(3, boundary::east,  condition_type::dirichlet, &g, 0);

    // --------------------------------------------------------------------------------

    /*
    // Homogeneous (zero) Dirichlet
    bcInfo.addCondition(0, boundary::north,  condition_type::dirichlet, &g_zero);
    bcInfo.addCondition(1, boundary::north,  condition_type::dirichlet, &g_zero);
    bcInfo.addCondition(2, boundary::north,  condition_type::dirichlet, &g_zero);
    bcInfo.addCondition(3, boundary::north,  condition_type::dirichlet, &g_zero);

    // Heterogeneous (non-zero) Dirichlet
    bcInfo.addCondition(0, boundary::south,  condition_type::dirichlet, &g_dir);
    bcInfo.addCondition(2, boundary::south,  condition_type::dirichlet, &g_dir);

    // Homogeneous (zero) Neumann
    bcInfo.addCondition(1, boundary::south, condition_type::neumann, &h_neu);
    bcInfo.addCondition(3, boundary::south, condition_type::neumann, &h_neu);
    */


    // Geometry case: Square
/*
    // ---------------Define Geometry---------------
    // (Unit square with 4 patches)
    //gsMultiPatch<> * patches;// = new gsMultiPatch<>;
    //patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);


    // ---------------Define Boundary conditions---------------
    //Dirichlet BCs
    bcInfo.addCondition(0, boundary::south,  condition_type::dirichlet, &g);
    bcInfo.addCondition(2, boundary::south,  condition_type::dirichlet, &g);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, &g);
    bcInfo.addCondition(2, boundary::east,  condition_type::dirichlet, &g);
    bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, &g);
    bcInfo.addCondition(3, boundary::east,  condition_type::dirichlet, &g);

    // Neumann BCs
    //h(x) = grad(u)*n + pn
    gsFunctionExpr<> h_neu("-2*pi*(1+0.25*cos(2*pi*y))","0.0", 2);
    bcInfo.addCondition(0, boundary::west, condition_type::neumann, &h_neu);
    bcInfo.addCondition(1, boundary::west, condition_type::neumann, &h_neu);

    // --------------------------------------------------------------------------------
*/

    dirichlet::strategy Dstrategy = dirichlet::elimination;
    iFace::strategy     Istrategy = iFace::glue;

    gsStokesPde<real_t> stokesPde(*patches, bcInfo, &f);

    gsMultiBasis<> multibasis (*patches);

    for (int i = 0; i < numRefine; ++i)
    {
        multibasis.uniformRefine();
    }

    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int max_tmp = multibasis.minCwiseDegree();

        // Elevate all degrees uniformly
        max_tmp += numElevate;
        multibasis.setDegree(max_tmp);
    }


    std::vector<gsBasis<>*> tmpBasisVec;
    for (size_t i = 0; i < patches->nPatches(); ++i)
    {
        tmpBasisVec.push_back(const_cast<gsBasis<>*>(&multibasis.basis(i)));
    }

    // if you need more test functions
    // std::for_each(tmpBasisVec.begin(), tmpBasisVec.end(), uniformRefine);


    // constructs appropriate spaces, which should produce good results
    // spaces satisfy the conditions which are described in Andrea's thesis
    std::vector<std::vector<gsBasis<>*> > bases;
    std::vector<gsPhysicalSpace*> phySpace = constructTHSpaces(tmpBasisVec,
                                                               *patches,
                                                               &bases);

    // --------------------------------------------------------------------------------
    // ˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇ
    // printing of the degrees of freedom
    int numDofs = 0;
    for (std::size_t i = 0; i != bases.size(); i++)
    {
        int numCompDofs = 0;
        std::cout << "i = " << i << "\n----------------------------------\n\n";
        
        for (std::size_t j = 0; j != bases[i].size(); j++)
        {
            numCompDofs += bases[i][j]->size();
        }

        numDofs += numCompDofs;
        std::cout << "Num comp dofs: " << numCompDofs << std::endl;
    }
    std::cout << "Num DOFS: " << numDofs << std::endl;
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // --------------------------------------------------------------------------------
    
    gsRecipeAssemblerStokes assembler(stokesPde);
    assembler.setSpace(phySpace);
    assembler.setZeroAverage(false);
    assembler.assemble();

    gsSparseMatrix<> sys = assembler.getSystemMatrix();
    gsMatrix<> rhs = assembler.getSystemRhs();

    gsMatrix<> eli  = solve(assembler.getEliminatedMatrix(),assembler.getEliminatedRhs());
    rhs -= assembler.getRhsModMatrix()*eli;

    gsMatrix<> sol = solve(sys,rhs);

    gsMatrix<> veloc_coefs = assembler.reconstructSolution(0, sol, eli);
    gsMatrix<> press_coefs = assembler.reconstructSolution(1, sol, eli);

    // --------------------------------------------------------------------------------
    // ˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇ
    // printing the sizes of the coefficients
    
    std::cout << "solution size: " << sol.rows() << " x " << sol.cols() << "\n"
              << "velocity size: " << veloc_coefs.rows() << " x " << veloc_coefs.cols() << "\n"
              << "pressure size: " << press_coefs.rows() << " x " << press_coefs.cols() << "\n"
              << std::endl;

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // --------------------------------------------------------------------------------


    // ================================================================================
    // reconstruction of the solution

    gsBoxTopology topology = getTopology(*patches);

    gsMultiBasis<> veloc_basis = getVelocityBasis(bases, topology);
    gsMultiBasis<> press_basis = getPressureBasis(bases, topology);

    gsMultiPatch<> velocityMp = buildVelocitySolution(veloc_basis, veloc_coefs, phySpace);
    gsMultiPatch<> pressureMp = buildPressureSolution(press_basis, press_coefs, phySpace);

    gsField<> velocity(*patches, velocityMp);
    gsField<> pressure(*patches, pressureMp);

    gsWriteParaview(velocity, output + "Velocity");
    gsWriteParaview(pressure, output + "Pressure");
    const gsField<> exact(*patches, g, false);
    gsWriteParaview<>( exact, output + "Exact", 1000 );

    
    return 0;
}
    


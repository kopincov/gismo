/** @file tutorialAssembly.cpp

    @brief Tutorial on how to use the gsAssembler class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

# include <gismo.h>
#include <iomanip>

#include <gsPde/gsPoissonHeterogeneousPde.h>
//#include <gsAssembler/gsVisitorPoissonHeterogeneous.h> // Stiffness volume integrals and load vector
#include <gsAssembler/gsVisitorDg2.h>
//#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
//#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals
//#include <gsAssembler/gsVisitorNitsche2.h> // Nitsche boundary integrals
#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>

using namespace gismo;

void reservemem(gsPoissonHeterogeneousAssembler<real_t> & assembler)
{
    gsSparseEntries<real_t> se;

    //estimate the number of triplets
    real_t OV = assembler.options().getReal("bdO"); //backup
    assembler.options().setReal("bdO", 0);   //set to zero
    se.reserve(assembler.system().cols() * assembler.numColNz() ); //for each column 2*p+1
    assembler.options().setReal("bdO", OV);//use the backup
    
    const gsMultiBasis<real_t>& mb = assembler.multiBasis(0);
    gsMatrix<index_t> actives;
    for(size_t np=0; np<assembler.pde().domain().nPatches();++np)
    {
        const gsBasis<real_t> &basis = mb[np];

        // Initialize domain element iterator -- using unknown 0
        gsBasis<>::domainIter domIt = basis.makeDomainIterator();

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            basis.active_into(domIt->centerPoint(),actives);
            assembler.system().mapColIndices(actives,np,actives);

            for(index_t i=0; i<actives.rows();i++)
                for(index_t j=0;j<actives.rows();j++)
                    se.addSorted(i,j,0);
        }
    }

    assembler.system().rhs().setZero(assembler.system().matrix().cols(), assembler.pde().numRhs());
    assembler.system().matrix().setFrom(se);
}



void printStr(std::string t, const int& width)
{
    gsInfo << std::left << std::setw(width) << std::setfill(' ') << t;
}

template<typename T> void printField(T t, const int& width, int prec = 3)
{
    gsInfo << std::left << std::setw(width) << std::setfill(' ') <<std::setprecision(prec);
    gsInfo << t;
}






int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool test3D=false;
    bool elevate = false;
    bool noSpPattern = false;
    index_t max_ref = 2;
    index_t max_deg = 1;


    gsCmdLine cmd("Checking the assembling speed for Poisson");
    cmd.addSwitch("d", "dim", "test3D", test3D);
    cmd.addSwitch("elev","Use elevation or increase", elevate);
    cmd.addSwitch("noSp","no test with sparsity pattern",noSpPattern);
    cmd.addInt("r","refine","Number of refinements",max_ref);
    cmd.addInt("p","degree","maximal degree",max_deg);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    int dim;
    test3D ? dim =3:dim =2;
    gsInfo<<"Dimension of Domain: "<<dim<<"\n";
    //! [Parse command line]


    // Grab a pre-defined grid of patches
    gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);
    gsMultiPatch<> patches3D;
    if(test3D)
    {
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&patches.patch(i));

            patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb));
        }
        patches3D.computeTopology();

        patches = give(patches3D);
    }

    gsInfo << "The domain is: "<< patches <<"\n";

    //! [Boundary conditions]
    gsFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)+pi/10",dim);

    gsBoundaryConditions<> bcInfo;

    for (gsMultiPatch<>::const_biterator
         bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
    }
    //! [Boundary conditions]

    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",dim);
    gsPiecewiseFunction<real_t>  alpha;
    gsFunctionExpr<> a1("1",dim);
    gsFunctionExpr<> a2("1",dim);
    for(size_t np = 0; np<patches.nPatches();np++)
        if(np%2 == 0)
            alpha.addPiece(a1);
        else
            alpha.addPiece(a2);


    gsPoissonHeterogeneousPde<real_t> ppde(patches, bcInfo, f, alpha);
    //! [Poisson Pde]

    // Copy basis from the multi-patch ( one per patch)
    gsMultiBasis<> splinebasis_backup( patches );


    gsVector<real_t> memOverheads(7);
    memOverheads<< 0, 0 , 0.1, 0.2 ,0.25 , 0.3, 0.5;

    gsOptionList opt = gsAssembler<>::defaultOptions();
    opt.setInt("DirichletValues"  , dirichlet::l2Projection);
    opt.setInt("DirichletStrategy", dirichlet::elimination );
    opt.setInt("InterfaceStrategy", iFace    ::conforming  );
    opt.setReal("quA", 1.0);
    opt.setInt ("quB", 1  );
    opt.setReal("bdA", 2.0);
    opt.setInt ("bdB", 1  );
    gsInfo << "Assembler "<< opt;

    gsMatrix<> totalAssembling(max_deg,memOverheads.rows());
    gsMatrix<> assembling(max_deg,memOverheads.rows());
    gsMatrix<> reserving(max_deg,memOverheads.rows());

    gsStopwatch t1,t2;
    //The mOv ==0 -> use sparsity pattern


    for(index_t mOv =0; mOv<memOverheads.rows();++mOv)
    {
        if(mOv==0&&noSpPattern)
        {
            totalAssembling.col(0).setZero();
            assembling.col(0).setZero();
            reserving.col(0).setZero();

            continue;
        }
        opt.setReal("bdO", memOverheads(mOv) );
        for(int deg=0; deg < max_deg;++deg)
        {
            gsMultiBasis<> splinebasis(splinebasis_backup);

            // Elevate all degrees uniformly
            if(elevate)
                splinebasis.degreeElevate(deg);
            else
                splinebasis.degreeIncrease(deg);

            // h-refine each basis (4, one for each patch)
            for (int i = 0; i < max_ref; ++i)
                splinebasis.uniformRefine();
            //! [Refinement]


            gsPoissonHeterogeneousAssembler<real_t>* PA = new gsPoissonHeterogeneousAssembler<real_t>();
            PA->initialize(ppde, splinebasis, opt);

            //! [Fill matrix manually]
            t2.restart();
            t1.restart();
            // reserve memeory
            if(mOv==0)
                reservemem(*PA);
            else
                PA->system().reserve(splinebasis[0] ,opt, ppde.numRhs());
            reserving(deg,mOv)=t1.stop();

            //start assembling
            t1.restart();

            PA->computeDirichletDofs();

            gsAssembler<real_t>* A = PA;
            A->push<gsVisitorPoissonHeterogeneous<real_t> >();

            A->push<gsVisitorNeumann<real_t> >(ppde.bc().neumannSides() );

            if ( opt.getInt("InterfaceStrategy") == iFace::dg )
                PA->push<gsVisitorDg2<real_t> >(ppde.domain().interfaces());

            PA->finalize();

            assembling(deg,mOv)=t1.stop();
            totalAssembling(deg,mOv)=t2.stop();
            index_t max=0;
            for(index_t i=0; i<PA->system().matrix().cols();++i)
                if(max<PA->system().matrix().col(i).nonZeros())
                    max = PA->system().matrix().col(i).nonZeros();
            gsInfo<<"Max number of entries: "<<max <<"\n"<<"usage: "
                  <<(real_t)max / PA->numColNz()<<"\n";
            //! [Fill matrix manually]

            gsInfo << "Assembled a system (matrix and load vector) with "
                   << PA->numDofs() << " dofs.\n";
            delete PA;
        }
    }

    //print the results;

    //reservation times
    gsInfo<<"\n Time for reserving the memory\n";
    printStr("deg",6);
    gsInfo<<"|";
    printStr("spPat",10);
    for(int mOv =1; mOv<memOverheads.rows();++mOv)
        printField(memOverheads(mOv),10);
    gsInfo<<"\n";
    gsInfo<<"-----------------------------------------------------------------\n";
    for(int deg=0; deg < max_deg;++deg)
    {
        printField(deg+1,6,0);
        gsInfo<<"|";
        for(index_t mOv =0; mOv<memOverheads.rows();++mOv)
            printField(reserving(deg,mOv),10);
        gsInfo<<"\n";
    }

    //assembling times
    gsInfo<<"\n Time for assembling\n";
    printStr("deg",6);
    gsInfo<<"|";
    printStr("spPat",10);
    for(int mOv =1; mOv<memOverheads.rows();++mOv)
        printField(memOverheads(mOv),10);
    gsInfo<<"\n";
    gsInfo<<"-----------------------------------------------------------------\n";
    for(int deg=0; deg < max_deg;++deg)
    {
        printField(deg+1,6,0);
        gsInfo<<"|";
        for(index_t mOv =0; mOv<memOverheads.rows();++mOv)
            printField(assembling(deg,mOv),10);
        gsInfo<<"\n";
    }

    //total times times
    gsInfo<<"\n Time for reserving and assembling\n";
    printStr("deg",6);
    gsInfo<<"|";
    printStr("spPat",10);
    for(index_t mOv =1; mOv<memOverheads.rows();++mOv)
        printField(memOverheads(mOv),10);
    gsInfo<<"\n";
    gsInfo<<"-----------------------------------------------------------------\n";
    for(int deg=0; deg < max_deg;++deg)
    {
        printField(deg+1,6,0);
        gsInfo<<"|";
        for(index_t mOv =0; mOv<memOverheads.rows();++mOv)
            printField(totalAssembling(deg,mOv),10);
        gsInfo<<"\n";
    }

    return 0;
}


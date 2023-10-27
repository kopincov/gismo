/**  gsParallelOperatorTest.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):      C. Hofer
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI

#include <gsIETI/gsParallelOperator.h>
#include <gismo.h>
#include <gsAssembler/gsStokesAssemblerNew.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsBoundary.h>

#include <gsMultiGrid/gsParallelGridHierarchy.h>
#include <gsMultiGrid/gsParallelMultiPatchPreconditioners.h>
#include <gsMultiGrid/gsMultiPatchPreconditioners.h>
#include <gsMultiGrid/gsMassSmoother.h>
#include <gsAssembler/gsParameterDomainAssembler.h>



using namespace gismo;

int main (int argc, char** argv)
{

    gsStopwatch time;
    bool pV = false;
    bool pVA = false;
    bool noSAOp = false;
    bool speedTest = false;
    int numRefine = 2;
    int nPatches  =  4;
    gsCmdLine cmd("An example for testing MPI with G+Smo.\n");
    cmd.addSwitch("",  "printV", "prints vectors", pV);
    cmd.addSwitch("",  "printVA", "prints all vectors", pVA);
    cmd.addInt   ("r", "uniformRefine",      "Number of uniform h-refinement steps to perform before testing",  numRefine        );
    cmd.addInt   ("n",  "nPatches",              "nPatches in each direction",       nPatches         );
    cmd.addSwitch("","NoSubassembledOperators", "does not use the subassembledOperator", noSAOp);
    cmd.addSwitch("", "SpeedTest","Perform several application of the operator to tests its serial performance",speedTest);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsOptionList options = cmd.getOptionList();
    
    if(pVA==true) pV=true;

    int nIt = 1;
    if(speedTest) nIt = 2000;

    // Initialize the MPI environment
    const gsMpi & mpi = gsMpi::init(argc, argv);

    // Get the world communicator
    gsMpiComm comm = mpi.worldComm();

    //Get size and rank of the processor
    int _size = comm.size();
    int _rank = comm.rank();

    if (0==_rank)
        gsInfo<<"Running on "<<_size<<" processes.\n";
    comm.barrier();

    gsInfo <<"MPI is "<< (mpi.initialized() ? "" : "NOT ")
          <<"initialized on process "<< _rank <<"\n";
    comm.barrier();

    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",2);
    gsFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)+pi/10",2);

    // Print out source function and solution
    if(_rank==0)gsInfo<<"Source function "<< f << "\n";
    if(_rank==0)gsInfo<<"Exact solution "<< g <<"\n\n";

    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;
    patches = gsNurbsCreator<>::BSplineSquareGrid(nPatches, nPatches, 0.5);

    gsBoundaryConditions<> bcInfo;
    for (gsMultiPatch<>::const_biterator  bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
        bcInfo.addCondition( *bit, condition_type::dirichlet, &g );


    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );

    int numElevate = 1;

    // h-refine each basis (4, one for each patch)
    for (int i = 0; i < numRefine; ++i)
        refine_bases.uniformRefine();

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int max_tmp = refine_bases.minCwiseDegree();
        max_tmp += numElevate;
        refine_bases.setDegree(max_tmp);
    }
    int degree = refine_bases.degree();

    unsigned ratio = (unsigned)(patches.nPatches()/comm.size());
    gsSortedVector<size_t> myPatches;
    for(size_t np =0; np< ratio*comm.size();++np)
        if(np /ratio == (unsigned)_rank)
            myPatches.push_sorted_unique(np);
    for(size_t np=ratio*comm.size(); np<patches.nPatches();++np)
        if(np % comm.size() == (unsigned)_rank)
            myPatches.push_sorted_unique(np);
    
    gsInfo<<"My Patches: ";
    for(size_t i=0; i<myPatches.size();++i )
        gsInfo<<" "<<myPatches[i];
    gsInfo<<"\n";
    
    comm.barrier();

    gsPoissonAssembler<real_t> assembler(patches,refine_bases,bcInfo,f,dirichlet::elimination, iFace::glue);
    std::vector<gsMultiBasis<real_t> > bases(1,refine_bases);
    gsSparseMatrix<real_t> mat;
    if(_rank==0)
    {
        assembler.assemble();
        mat = assembler.matrix();
    }
    /*
    gsStokesPde<real_t> pde(*patches, bcInfo,&f);
    std::vector<gsMultiBasis<real_t> > bases(3);
    bases[2] = refine_bases;
    refine_bases.degreeIncrease();
    bases[0] = refine_bases; bases[1] = refine_bases;
    gsStokesAssemblerNew<real_t> assembler( pde,bases,dirichlet::elimination, iFace::glue);
*/
    gsDofMapper map = assembler.system().colMapper(0);

    /*
    std::vector<gsDofMapper> locMappers(patches->nPatches());
    //for(size_t npi = 0; npi<myPatches.size();++npi)
    for(index_t np = 0; np<patches->nPatches();++np)
    {
        //    size_t np = myPatches[npi];
        gsMultiBasis<real_t> mb(refine_bases.basis(np));
        gsBoundaryConditions<real_t> bc;
        bcInfo.getConditionsForPatch(np,bc);
        locMappers[np] = mb.getMapper(dirichlet::elimination,iFace::none,bc,0,true);
    }
*/
    // std::vector<std::vector<gsDofMapper> > locMappers(1,std::vector<gsDofMapper>(patches->nPatches()));
    std::vector<std::vector<gsDofMapper> > locMappers(1,std::vector<gsDofMapper>(patches.nPatches()));
    //for(size_t npi = 0; npi<myPatches.size();++npi)
    for(int c=0; c<1;++c)
        for(size_t np = 0; np<patches.nPatches();++np)
        {
            //    size_t np = myPatches[npi];
            gsMultiBasis<real_t> mb(bases[c].basis(np));
            gsBoundaryConditions<real_t> bc;
            bcInfo.getConditionsForPatch(np,bc);
            locMappers[c][np] = mb.getMapper(dirichlet::elimination,iFace::none,bc,0,true);
            //locMappers[c][np] = c==2 ? mb.getMapper(dirichlet::none,iFace::none,bc,0,true) : mb.getMapper(dirichlet::elimination,iFace::none,bc,0,true) ;
        }
    typedef gsParallelOperator<real_t>::connectionPair connectionPair;
    typedef gsParallelOperator<real_t>::procDof procDof;

    std::vector<gsMatrixOp<gsSparseMatrix<real_t> >::Ptr > _localOps(patches.nPatches());
    std::vector<gsLinearOperator<real_t>::Ptr> localOps(patches.nPatches());
    for(size_t i=0; i<myPatches.size();++i)
    {
        gsBoundaryConditions<real_t> bcLoc;
        bcInfo.getConditionsForPatch(myPatches[i],bcLoc);
        gsPoissonAssembler<real_t> assemblerLoc(patches.patch(myPatches[i]),refine_bases.basis(myPatches[i]),bcLoc,f,dirichlet::elimination, iFace::glue);
        /*
        gsStokesPde<real_t>* pdeL = dynamic_cast<gsStokesPde<real_t>*>(pde.restrictToPatch(myPatches[i]));
        std::vector<gsMultiBasis<real_t> > basesL(3);
        basesL[2] = bases[2].basis(myPatches[i]); basesL[0] = bases[0].basis(myPatches[i]); basesL[1] = bases[1].basis(myPatches[i]);
        gsStokesAssemblerNew<real_t> assemblerLoc( *pdeL,basesL,dirichlet::elimination, iFace::glue);
        */
        assemblerLoc.assemble();
        // gsInfo<<"Local mat: \n"<<i<<": "<<assemblerLoc.system().matrix().toDense()<<"\n";
        _localOps[myPatches[i]] =  makeMatrixOp(assemblerLoc.system().matrix().moveToPtr());
        localOps[myPatches[i]] = _localOps[myPatches[i]]; //same lifetime
        //delete pdeL;
    }
    // gsMatrix<bool> bInteraction(3,3);
    // bInteraction(0,0) = true;bInteraction(1,1) = true;bInteraction(0,2) = true;bInteraction(2,0) = true;bInteraction(1,2) = true;bInteraction(2,1) = true;bInteraction(2,2) = true;
    gsMatrix<bool> bInteraction(1,1);
    bInteraction(0,0) = true;


    std::vector<gsDofMapper> vecMap(1,map);
    gsPatchSubassembledTopology<real_t>::Ptr top = gsPatchSubassembledTopology<real_t>::make(myPatches,bases,bInteraction,locMappers,vecMap);
    top->reorderLike(map);
    gsPatchSubassambledLocalOperator<real_t>::Ptr locOperator = gsPatchSubassambledLocalOperator<real_t>::make(localOps,top);

    std::vector<gsDofMapper> procLocalMapper(comm.size());
    gsDofMapper procLocalPatchMapper = locOperator->getTopologyRow()->getLocalPatchMapper();
    procLocalMapper[_rank] = locOperator->getTopologyRow()->getLocalMapper();
    gsPatchInterfaceConnections<real_t> patchConnection(myPatches,bases,bInteraction,locMappers,procLocalPatchMapper,map,comm);
    patchConnection.init();
    gsDofMapper procGlobMapper = patchConnection.generateProcGlobalMapper();
    if(_rank == 0) gsInfo<<"Creating Parallel Operator"<<"\n";
    gsParallelOperator<real_t>::Ptr parOp = gsParallelOperator<real_t>::make(patchConnection.getConnectionPairs(),locOperator,comm);
    if(_rank == 0) gsInfo<<"Creating Parallel Operator finished"<<"\n";
    comm.barrier();
    gsParallelGlobalLocalHandler handl(procGlobMapper,procLocalMapper,comm);
    comm.barrier();

    {
        //Test the usual operation
        gsMatrix<real_t> inp, globInp,out;
        gsVector<real_t> lin = gsVector<real_t>::LinSpaced(assembler.numDofs(),0,assembler.numDofs()-1);
        globInp.setZero(lin.rows(),1);
        globInp.col(0) = lin;
        handl.extractLocalVector(globInp,inp);
        comm.barrier();
        if(pVA)gsInfo<<"rank "<<_rank<<": local input vector: "<<inp.transpose()<<"\n";
        parOp->accumulateAfterApplication(true);
        if(_rank == 0) gsInfo<<"Apply Parallel Operator"<<"\n";
        comm.barrier();
        time.restart();

        for(int n =0; n<nIt;++n)
            parOp->apply(inp,out);
        if(_rank == 0 && speedTest) gsInfo<< "Time needed to complete apply: "<<time.stop()<<"\n";
        parOp->getConnectionHandlerTestSpace()->distributeAccumulatedVector(out);


        if(_rank == 0) gsInfo<<"Building global vector"<<"\n";
        gsMatrix<real_t> globVector;
        handl.buildGlobalVector(out,globVector);
        for(int i=0; i<comm.size();++i)
        {
            comm.barrier();
            if(_rank == 0)
            {
                if(pVA)gsInfo<<"Result for rank "<<_rank<<"\n"<<out.transpose()<<"\n\n"<<std::flush;
                if(pV)gsInfo<<"Full Result for rank "<<_rank<<"\n"<<globVector.transpose()<<"\n\n"<<std::flush;
            }
        }


        // Generate system matrix and load vector
        if(_rank==0)
        {
            if(pV)gsInfo<<"Total Result: \n"<<(mat*globInp).transpose()<<"\n"<<std::flush;
            gsInfo<<"|| A*u_MPI - A*u|| = "<<(mat*globInp-globVector).norm()<<"\n"<<std::flush;
            //gsInfo<< "Global Mat: \n"<<mat.toDense()<<"\n";


            gsMatrix<real_t> m;
            //locOperator->toMatrix(m);
            //gsInfo<< "Subass Mat: \n"<<m<<"\n";
        }
    }
    comm.barrier();
    if(_rank == 0) gsInfo<<"\n\n ------------     Start reduced apply  ------------"<<"\n\n"<<std::flush;

    gsBoundaryConditions<real_t> reducedBc;
    gsMultiPatch<real_t> reducedPatch;
    gsMultiBasis<real_t> reducedBasis;
    for(size_t i=0; i<myPatches.size();++i)
    {
        gsBoundaryConditions<real_t> bcLoc;
        bcInfo.getConditionsForPatch(myPatches[i],bcLoc);
        gsBoundaryConditions<real_t>::bcContainer container = bcLoc.allConditions();
        for(gsBoundaryConditions<real_t>::const_iterator it = container.begin(); it!=container.end();++it)
            reducedBc.add(i,it->side(),it->ctype(),it->function(),it->unknown(),it->unkComponent(),it->parametric());

        reducedPatch.addPatch(patches.patch(myPatches[i]).clone());
        reducedBasis.addBasis(bases.front().basis(myPatches[i]).clone().release());
    }
    reducedPatch.computeTopology();
    reducedBasis.setTopology(reducedPatch.topology());
    gsPoissonAssembler<real_t> reducedAssembler(reducedPatch,reducedBasis,reducedBc,f,dirichlet::elimination, iFace::glue);
    /*
    gsStokesPde<real_t>* pdeL = dynamic_cast<gsStokesPde<real_t>*>(pde.restrictToPatch(myPatches[i]));
    std::vector<gsMultiBasis<real_t> > basesL(3);
    basesL[2] = bases[2].basis(myPatches[i]); basesL[0] = bases[0].basis(myPatches[i]); basesL[1] = bases[1].basis(myPatches[i]);
    gsStokesAssemblerNew<real_t> assemblerLoc( *pdeL,basesL,dirichlet::elimination, iFace::glue);
    */



    gsPatchSubassembledTopology<real_t>::Ptr topNew = gsPatchSubassembledTopology<real_t>::make(myPatches,bases,bInteraction,locMappers,vecMap);
    topNew->reorderLike(map);
    procLocalPatchMapper = topNew->getLocalPatchMapper();
    procLocalMapper[_rank] = topNew->getLocalMapper();
    gsPatchInterfaceConnections<real_t> patchConnectionNew(myPatches,bases,bInteraction,locMappers,procLocalPatchMapper,map,comm);
    patchConnectionNew.init();
    procGlobMapper = patchConnectionNew.generateProcGlobalMapper();
    comm.barrier();
    gsParallelGlobalLocalHandler handlNew(procGlobMapper,procLocalMapper,comm);

    topNew->reorder(map,reducedAssembler.system().colMapper(0));
    reducedAssembler.assemble();

    gsSparseMatrix<real_t>& reducedMat = reducedAssembler.system().matrix();
    gsMatrixOp<gsSparseMatrix<real_t> >::Ptr reducedMatOp = makeMatrixOp(reducedMat.moveToPtr());
    gsParallelOperator<real_t>::Ptr reducedParOp = gsParallelOperator<real_t>::make(patchConnectionNew.getConnectionPairs(),reducedMatOp,comm);
    comm.barrier();

    {
        //Test the usual operation
        gsMatrix<real_t> inp, globInp,out;
        gsVector<real_t> lin = gsVector<real_t>::LinSpaced(assembler.numDofs(),0,assembler.numDofs()-1);
        globInp.setZero(lin.rows(),1);
        globInp.col(0) = lin;
        handlNew.extractLocalVector(globInp,inp);
        comm.barrier();
        if(pVA)gsInfo<<"rank "<<_rank<<": local input vector: "<<inp.transpose()<<"\n";
        if(_rank == 0) gsInfo<<"Apply Parallel Operator"<<"\n";
        comm.barrier();
        time.restart();

        for(int n =0; n<nIt;++n)
            reducedParOp->apply(inp,out);
        if(_rank == 0 && speedTest) gsInfo<< "Time needed to complete reduced apply: "<<time.stop()<<"\n";

        if(_rank == 0) gsInfo<<"Building global vector"<<"\n";
        gsMatrix<real_t> globVector;
        handlNew.buildGlobalVector(out,globVector);
        for(int i=0; i<comm.size();++i)
        {
            comm.barrier();
            if(_rank == 0)
            {
                if(pVA)gsInfo<<"Result for rank "<<_rank<<"\n"<<out.transpose()<<"\n\n"<<std::flush;
                if(pV)gsInfo<<"Full Result for rank "<<_rank<<"\n"<<globVector.transpose()<<"\n\n"<<std::flush;
            }
        }


        // Generate system matrix and load vector
        if(_rank==0)
        {
            if(pV)gsInfo<<"Total Result: \n"<<(mat*globInp).transpose()<<"\n"<<std::flush;
            gsInfo<<"|| A*u_MPI - A*u|| = "<<(mat*globInp-globVector).norm()<<"\n"<<std::flush;
            //gsInfo<< "Global Mat: \n"<<mat.toDense()<<"\n";

            /*
            gsMatrix<real_t> m;
            locOperator->toMatrix(m);
            gsInfo<< "Subass Mat: \n"<<m<<"\n\n";
            gsInfo<< "Reduced Mat: \n"<<reducedMat.toDense()<<"\n";
            */
        }
    }

    comm.barrier();
    if(_rank == 0) gsInfo<<"\n\n ------------     Start Hierarchy ------------"<<"\n\n"<<std::flush;

    gsParallelGridHierarchy<real_t> GH = gsParallelGridHierarchy<real_t>::buildByCoarsening(refine_bases,bcInfo,options,myPatches,mpi.worldComm());
    comm.barrier();
    if(_rank == 0) gsInfo<<"Build Restriction and Prolongation ops"<<"\n"<<std::flush;
    std::pair<std::pair<std::vector<gsParallelOperator<real_t>::Ptr>,std::vector<gsParallelOperator<real_t>::Ptr> >,
            std::vector<gsParallelGlobalLocalHandler::Ptr> >
            ops = GH.getRestrictionAndProlongationOperators();
    comm.barrier();
    {
        if(_rank == 0) gsInfo<<"Start Test"<<"\n"<<std::flush;
        gsMatrix<real_t> inpR,inpP, globInpR, globInpP,outR,outP;
        gsMatrix<real_t> globVectorR,globVectorP;

        gsParallelGlobalLocalHandler handlF = *ops.second[ops.second.size()-1];
        gsParallelGlobalLocalHandler handlC = *ops.second[ops.second.size()-2];
        size_t nDofsC = handlC.globalSize();
        size_t nDofsF = handlF.globalSize();
        if(_rank == 0) gsInfo<<"got handler: "<<nDofsC<<" - "<<nDofsF<<"\n"<<std::flush;
        comm.barrier();

        gsVector<real_t> lin = gsVector<real_t>::LinSpaced(nDofsF,0,nDofsF-1);
        globInpR.setZero(lin.rows(),1);
        globInpR.col(0) = lin;
        handlF.extractLocalVector(globInpR,inpR);

        if(_rank == 0) gsInfo<<"got extracted restricton vector"<<"\n"<<std::flush;
        comm.barrier();

        lin = gsVector<real_t>::LinSpaced(nDofsC,0,nDofsC-1);
        globInpP.setZero(lin.rows(),1);
        globInpP.col(0) = lin;
        handlC.extractLocalVector(globInpP,inpP);

        comm.barrier();
        if(_rank == 0) gsInfo<<"Restriction: Apply Parallel Operator"<<"\n"<<std::flush;
        ops.first.first.back()->getConnectionHandlerTrialSpace()->distributeAccumulatedVector(inpR);
        time.restart();
        for(int n =0; n<nIt;++n)
            ops.first.first.back()->apply(inpR,outR);
        if(_rank == 0 && speedTest) gsInfo<< "Time needed to complete restriction apply: "<<time.stop()<<"\n";
        handlC.buildGlobalVector(outR,globVectorR);

        comm.barrier();
        if(_rank == 0) gsInfo<<"Prolongation: Apply Parallel Operator"<<"\n"<<std::flush;
        for(int n =0; n<nIt;++n)
            ops.first.second.back()->apply(inpP,outP);
        if(_rank == 0 && speedTest) gsInfo<< "Time needed to complete prolongation apply: "<<time.stop()<<"\n";

        ops.first.second.back()->getConnectionHandlerTestSpace()->distributeAccumulatedVector(outP);
        handlF.buildGlobalVector(outP,globVectorP);

        for(int i=0; i<comm.size();++i)
        {
            comm.barrier();
            if(_rank == 0)
            {
                if(pVA)gsInfo<<"Result for rank "<<_rank<<" Restriction:  \n"<<outR.transpose()<<"\n\n"<<std::flush;
                if(pV)gsInfo<<"Full Result for rank "<<_rank<<" Restriction: \n"<<globVectorR.transpose()<<"\n\n"<<std::flush;
            }
        }

        for(int i=0; i<comm.size();++i)
        {
            comm.barrier();
            if(_rank == 0)
            {
                if(pVA)gsInfo<<"Result for rank "<<_rank<<" Prolongation: \n"<<outP.transpose()<<"\n\n"<<std::flush;
                if(pV)gsInfo<<"Full Result for rank "<<_rank<<" Prolongation: \n"<<globVectorP.transpose()<<"\n\n"<<std::flush;
            }
        }

        // Generate system matrix and load vector
        if(_rank==0)
        {


            gsInfo<< "Assembling...\n";
            gsGridHierarchy<real_t> GH_ = gsGridHierarchy<real_t>::buildByCoarsening(refine_bases,bcInfo,options);
            const gsSparseMatrix<real_t> & _mat = GH_.getTransferMatrices().back();
            const gsSparseMatrix<real_t> & _matT = GH_.getTransferMatrices().back().transpose();
            if(pV)gsInfo<<"Total Result Restriction: \n"<<(_matT*globInpR).transpose()<<"\n"<<std::flush;
            if(pV)gsInfo<<"Total Result Prolongation: \n"<<(_mat*globInpP).transpose()<<"\n"<<std::flush;
            gsInfo<<"|| R*u_MPI - R*u|| = "<<(_matT*globInpR-globVectorR).norm()<<"\n"<<std::flush;
            gsInfo<<"|| P*u_MPI - P*u|| = "<<(_mat*globInpP-globVectorP).norm()<<"\n"<<std::flush;
            // gsInfo<<"|| R*u_MPI - R*u|| = "<<(mat*globInp-globVector).norm()<<"\n";
            /*
            gsInfo<< "Global Mat: \n"<<matT.toDense()<<"\n\n";

            gsMatrix<real_t> m;
            ops.first.first.back()->toMatrix(m);
            gsInfo<< "Subass Mat: \n"<<m<<"\n\n\n";

            gsInfo<< "Global Mat: \n"<<mat.toDense()<<"\n\n";
            ops.first.second.back()->toMatrix(m);
            gsInfo<< "Subass Mat: \n"<<m<<"\n";
            */
        }
    }
    if(_rank == 0) gsInfo<<"------------------------------------------ Galerkin Projection of system matrices ------------------\n\n";

    std::vector<gsMatrixOp<gsSparseMatrix<real_t> >::Ptr > c_ops;
    gsMatrixOp<gsSparseMatrix<real_t> >::Ptr c_op;
    std::vector<gsDofMapper> c_locMappers;
    gsDofMapper c_globMapper;
    std::pair<std::vector<typename gsParallelOperator<real_t>::Ptr >, std::vector<gsParallelGlobalLocalHandler::Ptr> > operators;
    if(!noSAOp)
        operators= GH.generateGalerkinProjection(_localOps,c_ops, c_locMappers,c_globMapper);
    else
        operators= GH.generateGalerkinProjection(reducedMatOp,c_op, c_locMappers,c_globMapper);
    gsGridHierarchy<real_t> GH_ = gsGridHierarchy<real_t>::buildByCoarsening(refine_bases,bcInfo,options);
    comm.barrier();
    {
        if(_rank == 0) gsInfo<<"Start Test"<<"\n"<<std::flush;

        gsSparseMatrix<real_t> _mat = mat;
        for(int l=(int)operators.first.size()-1;l>=0; --l)
        {
            gsMatrix<real_t> inp,out, globInp,globVector;
            const gsParallelOperator<real_t>& op_l = *operators.first[l];
            const gsParallelGlobalLocalHandler& handler_l = *operators.second[l];
            size_t nDofs = handler_l.globalSize();
            if(_rank == 0) gsInfo<<"got handler: "<<nDofs<<"\n"<<std::flush;
            comm.barrier();

            gsVector<real_t> lin = gsVector<real_t>::LinSpaced(nDofs,0,nDofs-1);
            globInp.setZero(lin.rows(),1);
            globInp.col(0) = lin;
            handler_l.extractLocalVector(globInp,inp);

            if(_rank == 0) gsInfo<<"got local vector"<<"\n"<<std::flush;
            comm.barrier();

            if(_rank == 0) gsInfo<<"Apply parallel operator on level "<<l<<"\n"<<std::flush;

            time.restart();
            for(int n =0; n<nIt;++n)
                op_l.apply(inp,out);
            if(_rank == 0 && speedTest) gsInfo<< "Time needed to complete op on level "<<l<<": "<<time.stop()<<"\n";

            handler_l.buildGlobalVector(out,globVector);


            for(int i=0; i<comm.size();++i)
            {
                comm.barrier();
                if(_rank ==i)
                {
                    if(pVA)gsInfo<<"Result for rank "<<_rank<<" Matrix on level "<<l<<":  \n"<<inp.transpose()<<"\n\n"<<std::flush;
                    if(pV)gsInfo<<"Full Result for rank "<<_rank<<" Matrix on level "<<l<<": \n"<<globVector.transpose()<<"\n\n"<<std::flush;
                }
            }

            if(_rank==0)
            {
                if(l != (int)operators.first.size()-1)
                {
                    _mat = (GH_.getTransferMatrices()[l].transpose()*_mat*GH_.getTransferMatrices()[l]).eval();
                    //     gsDebug<<"\n"<<GH_.getTransferMatrices()[l].toDense()<<"\n\n";
                }
                //   gsDebug<<"matrix: \n"<<_mat.toDense()<<"\n\n";
                if(pV)gsInfo<<"Total Result level "<<l<<": \n"<<(_mat*globInp).transpose()<<"\n"<<std::flush;
                gsInfo<<"|| K*u_MPI - K*u||  level "<<l<<"= "<<(globVector-_mat*globInp).norm()<<"\n"<<std::flush;
            }


        }


    }

    if(_rank == 0) gsInfo<<"\n\n ------------     Start Additive precond------------"<<"\n\n"<<std::flush;

    parOp->accumulateAfterApplication(false);
    gsParallelAdditivePreconditionerOp<real_t>::uPtr addSmooth = gsParallelAdditivePreconditionerOp<real_t>::make(
                parOp,
                setupPiecewisePreconditioner(
                    *parOp,
                    handl,
                    makeSubspaceCorrectedMassSmootherOperatorsDirichlet(refine_bases,1),
                    patches,
                    refine_bases,
                    bcInfo,
                    options,
                    myPatches,
                    comm
                    ),
                1
                );

    comm.barrier();
    {
        if(_rank == 0) gsInfo<<"Start Test"<<"\n"<<std::flush;
        gsMatrix<real_t> inp, globInp,globVector;

        size_t nDofs = handl.globalSize();
        if(_rank == 0) gsInfo<<"got handler: "<<nDofs<<"\n"<<std::flush;
        comm.barrier();

        gsVector<real_t> lin = gsVector<real_t>::LinSpaced(nDofs,0,nDofs-1);
        globInp.setZero(lin.rows(),1);
        globInp.col(0) = lin;
        handl.extractLocalVector(globInp,inp);

        if(_rank == 0) gsInfo<<"got local vector"<<"\n"<<std::flush;
        comm.barrier();

        if(_rank == 0) gsInfo<<"Apply Parallel Smoother"<<"\n";
        time.restart();
        gsMatrix<real_t> inp_;
        for(int n =0; n<nIt;++n)
        {
            inp_ = inp;
            addSmooth->step(gsMatrix<real_t>::Zero(inp.rows(),inp.cols()),inp_);
        }
        inp = inp_;
        if(_rank == 0 && speedTest) gsInfo<< "Time needed to complete smoother step: "<<time.stop()<<"\n";

        parOp->getConnectionHandlerTestSpace()->distributeAccumulatedVector(inp);
        handl.buildGlobalVector(inp,globVector);


        for(int i=0; i<comm.size();++i)
        {
            comm.barrier();
            if(_rank == i)
            {
                if(pVA)gsInfo<<"Result for rank "<<_rank<<" Smoother:  \n"<<inp.transpose()<<"\n\n"<<std::flush;
                if(pV)gsInfo<<"Full Result for rank "<<_rank<<" Smoother: \n"<<globVector.transpose()<<"\n\n"<<std::flush;
            }
        }

        // Generate system matrix and load vector
        if(_rank==0)
        {


            gsInfo<< "Assembling...\n";
            gsAdditivePreconditionerOp<real_t>::Ptr sm = gsAdditivePreconditionerOp<real_t>::make(makeMatrixOp(mat),
                                                                                                  setupPiecewisePreconditioner(
                                                                                                      mat,
                                                                                                      makeSubspaceCorrectedMassSmootherOperatorsDirichlet(refine_bases,1),
                                                                                                      patches,
                                                                                                      refine_bases,
                                                                                                      bcInfo,
                                                                                                      options
                                                                                                      ),
                                                                                                  1);
            sm->step(gsMatrix<real_t>::Zero(globInp.rows(),globInp.cols()),globInp);

            if(pV)gsInfo<<"Total Result Smoother: \n"<<globInp.transpose()<<"\n"<<std::flush;
            gsInfo<<"|| S*u_MPI - S*u|| = "<<(globInp-globVector).norm()<<"\n"<<std::flush;
            /*
            gsInfo<< "Global Mat: \n"<<matT.toDense()<<"\n\n";

            gsMatrix<real_t> m;
            ops.first.first.back()->toMatrix(m);
            gsInfo<< "Subass Mat: \n"<<m<<"\n\n\n";

            gsInfo<< "Global Mat: \n"<<mat.toDense()<<"\n\n";
            ops.first.second.back()->toMatrix(m);
            gsInfo<< "Subass Mat: \n"<<m<<"\n";
            */
        }
    }


    if(_rank == 0) gsInfo<<"\n\n ------------     Start Additive precond 2------------"<<"\n\n"<<std::flush;


    gsParallelAdditivePreconditionerOp<real_t>::Ptr addSmooth2;
    {
        real_t damping =1;
        const index_t dim = refine_bases.dim();
        gsMultiBasis<> proc_loc_mb;
        for (size_t k=0;k<myPatches.size();++k)
            proc_loc_mb.addBasis(refine_bases[myPatches[k]].clone().release());

        const index_t proc_loc_dofs = (*operators.second.back()).localSize();

        gsDofMapper dm; // Global mapper
        refine_bases.getMapper(dirichlet::elimination, iFace::glue,  bcInfo,   dm, 0);

        // Global => Processor local
        gsMatrix<unsigned> glob2loc;
        {
            gsMatrix<unsigned> inp(proc_loc_dofs,1);
            inp.col(0) = gsVector<unsigned>::LinSpaced(proc_loc_dofs,0,proc_loc_dofs-1);
            glob2loc.setZero((*operators.second.back()).globalSize(),1);
            (*operators.second.back()).addLocalVectorToGlobal( inp, glob2loc );
        }

        std::vector< gsVector<index_t> > proc_loc_maps;
        for (size_t k=0;k<myPatches.size();++k)
        {
            const index_t sz = proc_loc_mb[k].size();
            index_t kk = myPatches[k];
            gsVector<index_t> local(sz);
            for (index_t j=0; j<sz; ++j)
            {
                local[j] = dm.is_free(j,kk) ? glob2loc(dm.index(j, kk),0) : -1u;
            }
            proc_loc_maps.push_back(give(local));
        }
        std::vector< std::vector< std::pair< typename gsBasis<>::Ptr, gsSparseMatrix<> > > > pieces =
                constructPieces( proc_loc_mb, proc_loc_maps, proc_loc_dofs );

        std::vector< gsLinearOperator<>::Ptr > localSmoothers;
        std::vector< gsSparseMatrix<> > smootherTransfers;
        const index_t nrPieces = pieces.size();

        real_t h = 1;
        for ( size_t j=0; j<refine_bases.nBases(); ++j)
            h = std::min(h,refine_bases[j].getMinCellLength());

        for ( index_t dd = 0; dd<nrPieces; ++dd )
        {
            gsBoundaryConditions<> localbc;
            for( index_t ps=0; ps < 2*dim; ++ps )
                localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

            const index_t sz = pieces[dd].size();
            for ( index_t j=0; j<sz; ++j)
            {
                if (dd == dim)
                {
                    localSmoothers.push_back( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][j].first, damping, localbc) );
                }
                else if (dd > 1)
                {
                    const real_t scalingK = std::pow(2,dim-dd) * (dim-dd) * std::pow( h/(2*degree+1), dim-dd );
                    const real_t scalingM = std::pow(2,dim-dd) * std::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                    localSmoothers.push_back( gsScaledOp<>::make( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][j].first, damping, localbc, scalingM/scalingK ), 1/scalingK ) );
                }
                else if (dd == 1)
                {
                    //gsSparseMatrix<> localMat = pieces[dd][j].second * mg.matrix(i) * pieces[dd][j].second.transpose();
                    //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );

                    const real_t scalingK = std::pow(2,dim-dd) * (dim-dd) * std::pow( h/(2*degree+1), dim-dd );
                    const real_t scalingM = std::pow(2,dim-dd) * std::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                    gsSparseMatrix<> M, K;
                    assembleParameterMass(*pieces[dd][j].first, M);
                    assembleParameterStiffness(*pieces[dd][j].first, K);
                    M = M.block( 1, 1, M.rows()-2, M.cols()-2 );
                    K = K.block( 1, 1, K.rows()-2, K.cols()-2 );
                    gsSparseMatrix<> mat = scalingK * K + scalingM * M;

                    localSmoothers.push_back( makeSparseCholeskySolver(mat) );
                }
                else if (dd == 0)
                {
                    //gsSparseMatrix<> localMat = pieces[dd][j].second * mg.matrix(i) * pieces[dd][j].second.transpose();
                    //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );
                    const index_t sz = pieces[dd][j].second.rows();
                    const real_t scalingFactor = std::pow(2,dim) * dim * std::pow( h/(1+2*degree), dim-1 ) * degree*degree/(h*(2*degree-1));
                    localSmoothers.push_back( gsScaledOp<>::make( gsIdentityOp<>::make(sz), 1/scalingFactor) );
                }

                smootherTransfers.push_back( give( pieces[dd][j].second ) );
            }
        }
        addSmooth2=give(
                    gsParallelAdditivePreconditionerOp<real_t>::make(
                        operators.first.back(),
                        give(smootherTransfers),
                        give(localSmoothers),
                        1
                        ));
    }



    comm.barrier();
    {
        if(_rank == 0) gsInfo<<"Start Test"<<"\n"<<std::flush;
        gsMatrix<real_t> inp, globInp,globVector;

        size_t nDofs = handl.globalSize();
        if(_rank == 0) gsInfo<<"got handler: "<<nDofs<<"\n"<<std::flush;
        comm.barrier();

        gsVector<real_t> lin = gsVector<real_t>::LinSpaced(nDofs,0,nDofs-1);
        globInp.setZero(lin.rows(),1);
        globInp.col(0) = lin;
        handl.extractLocalVector(globInp,inp);

        if(_rank == 0) gsInfo<<"got local vector"<<"\n"<<std::flush;
        comm.barrier();

        if(_rank == 0) gsInfo<<"Apply Parallel Smoother"<<"\n";
        time.restart();
        gsMatrix<real_t> inp_;
        for(int n =0; n<nIt;++n)
        {
            inp_ = inp;
            addSmooth2->step(gsMatrix<real_t>::Zero(inp.rows(),inp.cols()),inp_);
        }
        inp = inp_;
        if(_rank == 0 && speedTest) gsInfo<< "Time needed to complete smoother step: "<<time.stop()<<"\n";

        parOp->getConnectionHandlerTestSpace()->distributeAccumulatedVector(inp);
        handl.buildGlobalVector(inp,globVector);


        for(int i=0; i<comm.size();++i)
        {
            comm.barrier();
            if(_rank == i)
            {
                if(pVA)gsInfo<<"Result for rank "<<_rank<<" Smoother:  \n"<<inp.transpose()<<"\n\n"<<std::flush;
                if(pV)gsInfo<<"Full Result for rank "<<_rank<<" Smoother: \n"<<globVector.transpose()<<"\n\n"<<std::flush;
            }
        }

        // Generate system matrix and load vector
        if(_rank==0)
        {

            real_t damping =1;
            gsInfo<< "Assembling...\n";
            const index_t dim = refine_bases.dim();

            std::vector< std::vector< std::pair< typename gsBasis<>::Ptr, gsSparseMatrix<> > > > pieces =
                    constructPieces( refine_bases, bcInfo, options );

            std::vector< gsLinearOperator<>::Ptr > localSmoothers;
            std::vector< gsSparseMatrix<> > smootherTransfers;
            const index_t nrPieces = pieces.size();


            real_t h = 1;
            for ( size_t j=0; j<refine_bases.nBases(); ++j)
                h = std::min(h,refine_bases[j].getMinCellLength());

            for ( index_t dd = 0; dd<nrPieces; ++dd )
            {
                gsBoundaryConditions<> localbc;
                for( index_t ps=0; ps < 2*dim; ++ps )
                    localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

                const index_t sz = pieces[dd].size();
                for ( index_t j=0; j<sz; ++j)
                {
                    if (dd == dim)
                    {
                        localSmoothers.push_back( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][j].first, damping, localbc) );
                    }
                    else if (dd > 1)
                    {
                        const real_t scalingK = std::pow(2,dim-dd) * (dim-dd) * std::pow( h/(2*degree+1), dim-dd );
                        const real_t scalingM = std::pow(2,dim-dd) * std::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                        localSmoothers.push_back( gsScaledOp<>::make( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][j].first, damping, localbc, scalingM/scalingK ), 1/scalingK ) );
                    }
                    else if (dd == 1)
                    {
                        //gsSparseMatrix<> localMat = pieces[dd][j].second * mg.matrix(i) * pieces[dd][j].second.transpose();
                        //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );

                        const real_t scalingK = std::pow(2,dim-dd) * (dim-dd) * std::pow( h/(2*degree+1), dim-dd );
                        const real_t scalingM = std::pow(2,dim-dd) * std::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                        gsSparseMatrix<> M, K;
                        assembleParameterMass(*pieces[dd][j].first, M);
                        assembleParameterStiffness(*pieces[dd][j].first, K);
                        M = M.block( 1, 1, M.rows()-2, M.cols()-2 );
                        K = K.block( 1, 1, K.rows()-2, K.cols()-2 );
                        gsSparseMatrix<> mat = scalingK * K + scalingM * M;

                        localSmoothers.push_back( makeSparseCholeskySolver(mat) );
                    }
                    else if (dd == 0)
                    {
                        //gsSparseMatrix<> localMat = pieces[dd][j].second * mg.matrix(i) * pieces[dd][j].second.transpose();
                        //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );
                        const index_t sz = pieces[dd][j].second.rows();
                        const real_t scalingFactor = std::pow(2,dim) * dim * std::pow( h/(1+2*degree), dim-1 ) * degree*degree/(h*(2*degree-1));
                        localSmoothers.push_back( gsScaledOp<>::make( gsIdentityOp<>::make(sz), 1/scalingFactor) );
                    }

                    smootherTransfers.push_back( give( pieces[dd][j].second ) );
                }
            }
            gsAdditiveSmoother::Ptr sm =
                    gsAdditiveSmoother::make(
                        makeMatrixOp(mat),
                        give(smootherTransfers),
                        give(localSmoothers),
                        1
                        );

            sm->step(gsMatrix<real_t>::Zero(globInp.rows(),globInp.cols()),globInp);

            if(pV)gsInfo<<"Total Result Smoother: \n"<<globInp.transpose()<<"\n"<<std::flush;
            gsInfo<<"|| S*u_MPI - S*u|| = "<<(globInp-globVector).norm()<<"\n"<<std::flush;
            /*
            gsInfo<< "Global Mat: \n"<<matT.toDense()<<"\n\n";

            gsMatrix<real_t> m;
            ops.first.first.back()->toMatrix(m);
            gsInfo<< "Subass Mat: \n"<<m<<"\n\n\n";

            gsInfo<< "Global Mat: \n"<<mat.toDense()<<"\n\n";
            ops.first.second.back()->toMatrix(m);
            gsInfo<< "Subass Mat: \n"<<m<<"\n";
            */
        }
    }







    if(_rank == 0) gsInfo<<"------------------------------------------ Coarse Solver ------------------\n\n"<<std::flush;

    gsCoarseSolverAdapter<real_t>::Ptr cs;
    if(!noSAOp)
        cs= give(constructCoarseSolver<real_t>(c_ops,operators.second.front(),myPatches,c_locMappers,c_globMapper));
    else
        cs = give(constructCoarseSolver<real_t>(c_op,operators.second.front()));
    comm.barrier();
    {
        gsSparseMatrix<real_t> _mat = mat;
        if(_rank == 0) gsInfo<<"Start Test"<<"\n"<<std::flush;

        gsMatrix<real_t> inp,inp_,out,out_, globInp,globVector;
        const gsParallelOperator<real_t>& op_l = *operators.first.front();
        const gsParallelGlobalLocalHandler& handler_l = *operators.second.front();
        size_t nDofs = handler_l.globalSize();
        if(_rank == 0) gsInfo<<"got handler: "<<nDofs<<"\n"<<std::flush;
        comm.barrier();

        gsVector<real_t> lin = gsVector<real_t>::LinSpaced(nDofs,0,nDofs-1);
        globInp.setZero(lin.rows(),1);
        globInp.col(0) = lin;
        handler_l.extractLocalVector(globInp,inp_);
        op_l.distribute(inp_,inp); // the input needs to be distributed

        if(_rank == 0) gsInfo<<"got local vector"<<"\n"<<std::flush;
        comm.barrier();

        if(_rank == 0) gsInfo<<"Apply parallel Coarse grid solver"<<"\n"<<std::flush;

        time.restart();
        for(int n =0; n<nIt;++n)
            cs->apply(inp,out_);
        if(_rank == 0 && speedTest) gsInfo<< "Time needed to complete coarse grid solver: "<<time.stop()<<"\n";

        op_l.distribute(out_,out);// the output is accumulated
        handler_l.buildGlobalVector(out,globVector);


        for(int i=0; i<comm.size();++i)
        {
            comm.barrier();
            if(_rank ==i)
            {
                if(pVA)gsInfo<<"Result for rank "<<_rank<<" Matrix" <<":  \n"<<inp.transpose()<<"\n\n"<<std::flush;
                if(pV)gsInfo<<"Full Result for rank "<<_rank<<" Matrix" <<": \n"<<globVector.transpose()<<"\n\n"<<std::flush;
            }
        }

        if(_rank==0)
        {
            for(int l=(int)operators.first.size()-2;l>=0; --l)
            {
                _mat = (GH_.getTransferMatrices()[l].transpose()*_mat*GH_.getTransferMatrices()[l]).eval();
            }
            gsLinearOperator<real_t>::Ptr solv = makeSparseCholeskySolver(_mat);
            gsMatrix<real_t> output;
            solv->apply(globInp,output);
            //     gsDebug<<"\n"<<GH_.getTransferMatrices()[l].toDense()<<"\n\n";

            //   gsDebug<<"matrix: \n"<<_mat.toDense()<<"\n\n";
            if(pV)gsInfo<<"Total Result"<<": \n"<<output.transpose()<<"\n"<<std::flush;
            gsInfo<<"|| CS*u_MPI - CS*u||  = "<<(globVector-output).norm()<<"\n"<<std::flush;
        }
    }


    return 0;
}

#else
#include <gismo.h>
using namespace gismo;
int main(int argc, char **argv)
{
    gsInfo<<" No Mpi enabled! \n";
    return 0;
}

#endif

/**  gsSpaceTimeSlicesMG.cpp

    This file is a test of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):     C. Hofer
    Created on:  2017-07-27

    DESCRIBE THE PURPOSE OF THE EXAMPLE HERE

    WHAT IS HAPPENING

    EXPECTED RESULT
*/

#include <gsCore/gsConfig.h>
#ifdef GISMO_WITH_MPI

#include <gismo.h>
#include <gismo_dev.h>

#include <gsMpi/gsMpi.h>

#include <gsAssembler/gsMGSpaceTimeAdapter.h>
#include <gsAssembler/gsSpaceTimesliceAssembler.h>
#include <gsPde/gsSpaceTimePoissonPde.h>
#include <gsAssembler/gsPoissonHeterogeneousAssembler.h>
#include <gsAssembler/gsSpaceTimeNorm.h>
#include <gsMultiGrid/gsGridHierarchy.h>
#include <gsSolver/gsTimeParallelMultigrid.h>
#include <gsSolver/gsConjugateGradient.h>
#include <gsCore/gsConnectedBoundaryComponent.h>


using namespace gismo;

void printNonZeros(gsSparseMatrix<real_t>& mat, std::string filename)
{
    std::fstream out(filename,  std::ofstream::out | std::ofstream::trunc);
    for (int k=0; k<mat.outerSize(); ++k)
        for (gsSparseMatrix<real_t>::InnerIterator it(mat,k); it; ++it)
        {
            if(math::abs(it.value())>1.e-10)
                out<<util::to_string(it.col())+" "+util::to_string(-it.row())+" "+util::to_string(it.value())+"\n";
        }
    out.close();
}

void doExperiment(const gsSpaceTimePoissonPde<real_t>& pde,gsMultiBasis<real_t> refine_bases, const gsOptionList& opt, gsFunction<real_t>& exact)
{
    gsStopwatch time;
    gsMatrix<real_t> solVMG;
    gsHeatSTSlapOperator<real_t>::HeatSTSlapSettings HeatSettings(opt.getGroup("TMGSlap"));

    bool detailedPrint = opt.getSwitch("detailedOutput");
    bool useIETI = opt.getInt("TMGSlap.DecompositionSolver") == gsHeatSTSlapOperator<real_t>::IETI;

    // h-refine each basis (4, one for each patch)
    for(int d=0; d<pde.domain().parDim()-1;++d)
        for (int i = 0; i < opt.getInt("refX"); ++i)
            refine_bases.uniformRefineComponent(d,1,1);

    for (int i = 0; i < opt.getInt("refT"); ++i)
        refine_bases.uniformRefineComponent(pde.domain().parDim()-1,1,1);

  //  typedef gsSpaceTimesliceAssembler<real_t>::Permutation Permutation;
    typedef gsSpaceTimesliceAssembler<real_t>::TPData TPData;
    gsSpaceTimesliceAssembler<real_t> assembler(pde,refine_bases,opt.getReal("theta"), dirichlet::elimination);

    /*
    {
    gsMultiBasis<real_t> tempBasis = refine_bases;

    gsMatrix<real_t> condsX(7,7+1);
    gsMatrix<real_t> condsY(7,7+1);
    for(int r=0; r<condsX.rows();++r)
    {
        gsMultiBasis<real_t> tempBasis2 = tempBasis;
        for(int e=1; e<condsX.cols();++e)
        {
            tempBasis2.degreeIncrease();
            gsSpaceTimesliceAssembler<real_t> SP_ASS(pde,tempBasis2,opt.getReal("theta"), dirichlet::elimination);
            SP_ASS.assemble(true,false);

            gsSpaceTimesliceAssembler<real_t>::TPData data = SP_ASS.getTensorProductData();
            Eigen::GeneralizedEigenSolver<Eigen::Matrix<real_t,Dynamic,Dynamic> > solv;
            solv.compute(data.Kt[0]->toDense(),data.Mt[0]->toDense(),true);

            Eigen::JacobiSVD<Eigen::Matrix<std::complex<real_t>,Dynamic,Dynamic> > svdX(solv.eigenvectors());
            condsX(r,e) = svdX.singularValues()(0)/ svdX.singularValues()(svdX.singularValues().size()-1);

            Eigen::JacobiSVD<Eigen::Matrix<std::complex<real_t>,Dynamic,Dynamic> > svdY(data.Mt[0]->toDense()*solv.eigenvectors());
            condsY(r,e) = svdY.singularValues()(0)/ svdY.singularValues()(svdY.singularValues().size()-1);
        }
        condsX(r,0) = tempBasis.basis(0).component(pde.domain().parDim()-1).size();
        condsY(r,0) = tempBasis.basis(0).component(pde.domain().parDim()-1).size();
        tempBasis.uniformRefineComponent(pde.domain().parDim()-1,1,1);

    }
    gsInfo<<"Condition numbers of X: \n"<<condsX<<"\n\n"<<std::flush;
    gsInfo<<"Condition numbers of Y: \n"<<condsY<<"\n\n"<<std::flush;
    return 0;
    }
*/
    //gsSparseMatrix<real_t>::Ptr STmat_fine; //= memory::make_shared(new gsSparseMatrix<real_t>(assembler.matrix().twistedBy(per))); //bring in TP form

    int nTimeSlices = assembler.getNTimeSlices();
    int timeLevels = (int)(log2(nTimeSlices))+1;
    int size  = gsMpi::worldSize();
    int rank = gsMpi::worldRank();
    int nCoresTime = opt.getInt("TimeCores");
    if(opt.getInt("SolverType") != 0 && opt.getInt("TimeCores") != size)
    {
        nCoresTime = size;
        if(rank == 0)
            std::cout << "Set the number of time cores to " << size << std::endl;
    }
    int nCoresSpace = 1;

    if(rank==0)
    {
        std::cout << "MPI threads:    " << size << std::endl;
        if(getenv("OMP_NUM_THREADS")!=NULL) std::cout << "OpenMP threads: " << atoi(getenv("OMP_NUM_THREADS")) << std::endl;
        else std::cout << "warning: OMP_NUM_THREADS not set!" << std::endl;
    }

    if(size<nCoresSpace*nCoresTime || nTimeSlices < nCoresTime)
    {
        if(rank==0) std::cout << "The number of cores does not fit!" << std::endl;
        return;
    }
    int nCores = nCoresSpace*nCoresTime;


    MPI_Group global_group;
    MPI_Comm_group(gsMpi::worldComm(), &global_group);
    MPI_Comm global_comm = gsMpi::worldComm();
    int ranksAll[nCores]; MPI_Group all_group; MPI_Comm all_comm;
    for(int i=0; i<nCores; i++) ranksAll[i] = i;
    MPI_Group_incl(global_group, nCores, ranksAll, &all_group);
    MPI_Comm_create(global_comm, all_group, &all_comm);

    if(rank<nCores)
    {

        //Setup MG to get the distribution of time processors
        gsGridHierarchy<real_t> gridH =  gsGridHierarchy<real_t>::buildByCoarsening(assembler.getSpaceBases(0),assembler.getSpaceBC(),opt.getGroup("MG"));
        //if having IETI
        const std::vector<gsMultiBasis<real_t> >& bases = gridH.getMultiBases();
        int nPatches = bases[0].nBases();
        std::vector<std::vector<gsSparseMatrix<real_t,RowMajor> > > PatchRestr(nPatches);
        for(int np=0; np<nPatches;++np)
        {
            gsBoundaryConditions<real_t> bc;
            assembler.getSpaceBC().getConditionsForPatch(np,bc);
            PatchRestr[np] = gsGridHierarchy<real_t>::buildByCoarsening(assembler.getSpaceBases(0).basis(np),bc,opt.getGroup("MG")).getTransferMatrices();
        }
        std::vector<gsSparseMatrix<real_t,RowMajor> > restr = gridH.getTransferMatrices();
        std::vector<gsSparseMatrix<real_t,RowMajor>::Ptr > prolongPtr(restr.size());
        for(size_t i=0; i<restr.size();++i)
        {
            prolongPtr[i] = restr[i].moveToPtr();
        }
        int spaceLevels = (int)restr.size()+1;
        STMultigridMPI<gsHeatSTSlapOperator<real_t> >::MPI_Comm_Vars MPICommVarsTime; MPICommVarsTime.size = size; MPICommVarsTime.my_rank = rank; MPICommVarsTime.Comm = all_comm;
        //    MPI_Barrier(time_comm);
        //    gsInfo<<"Rank "<<rank <<" finished transfermats.\n"<<std::flush;
        //Basic Setup of ST MG
        STMultigridMPI<gsHeatSTSlapOperator<real_t> > testMG(MPICommVarsTime, timeLevels, spaceLevels, opt.getGroup("TMG"));
        //    MPI_Barrier(time_comm);
        //    gsInfo<<"Rank "<<rank <<" finished STMG setup.\n"<<std::flush;

        //Do the space assembling
        time.restart();
        assembler.assemble(true,true);
        real_t assTime = time.stop();
        MPI_Barrier(all_comm);
        if(rank==0) assembler.printTimings();
        //  MPI_Barrier(time_comm);
        //  gsInfo<<"Rank "<<rank <<" finished assembling.\n"<<std::flush;

        //Setup of SlapOperator for deriving the coarsening patch
        std::vector<gsHeatSTSlapOperator<real_t>::Ptr > STOptemp(spaceLevels);
        if(rank==0)gsInfo<<"tau: "<<assembler.getTau(0)<<"\n";

        HeatSettings.m_spatialAssemblers.resize(spaceLevels);
        HeatSettings.m_STAssemblers.resize(spaceLevels);
        HeatSettings.m_pde.resize(spaceLevels);
        gsSpaceTimePoissonPde<real_t>::Ptr SlicePde = memory::make_shared(new gsSpaceTimePoissonPde<real_t>(*assembler.getSlice(0),assembler.getSpaceBC(),*pde.rhs(),*pde.rhs_x(),*pde.rhs_t(),*pde.getAlpha()));
        for(int step=0; step<spaceLevels; step++)
        {
            real_t h = assembler.getH(0)*math::exp2(spaceLevels-step-1);
            if(rank==0)gsInfo<<"h on level "<<step<<" :" <<h<<"\n";
            STOptemp[step] = memory::make_shared(new gsHeatSTSlapOperator<real_t>(h, assembler.getTau(0), timeLevels));
            testMG.setSTSlapOperator(STOptemp[step], step);
            HeatSettings.m_pde[spaceLevels-step-1] = assembler.getSpatialPde();
            HeatSettings.m_spatialAssemblers[spaceLevels-step-1] = memory::make_shared(new gsPoissonHeterogeneousAssembler<real_t>(*dynamic_cast<gsPoissonHeterogeneousPde<real_t>*>(HeatSettings.m_pde[spaceLevels-step-1].get()),bases[spaceLevels-step-1],dirichlet::elimination,iFace::glue));

            gsMultiBasis<real_t>::uPtr sliceBasis = memory::make_unique(new gsMultiBasis<real_t>());
            for(size_t i = 0; i<bases[spaceLevels-step-1].nBases();++i)
            {
                const gsBasis<real_t>& basis_i = bases[spaceLevels-step-1].basis(i);
                std::vector<gsBasis<real_t>*> basesVec;
                if(bases[spaceLevels-step-1].basis(i).domainDim()==1)
                {
                    basesVec.push_back(basis_i.component(0).clone().release());
                    basesVec.push_back(assembler.getTimeBases().basis(0).clone().release());
                    sliceBasis->addBasis(gsTensorBSplineBasis<2,real_t>::make(basesVec));
                }
                else if(bases[spaceLevels-step-1].basis(i).domainDim()==2)
                {
                    basesVec.push_back(basis_i.component(0).clone().release());
                    basesVec.push_back(basis_i.component(1).clone().release());
                    basesVec.push_back(assembler.getTimeBases().basis(0).clone().release());
                    sliceBasis->addBasis(gsTensorBSplineBasis<3,real_t>::make(basesVec));
                }
                else if(bases[spaceLevels-step-1].basis(i).domainDim()==3)
                {
                    basesVec.push_back(basis_i.component(0).clone().release());
                    basesVec.push_back(basis_i.component(1).clone().release());
                    basesVec.push_back(basis_i.component(2).clone().release());
                    basesVec.push_back(assembler.getTimeBases().basis(0).clone().release());
                    sliceBasis->addBasis(gsTensorBSplineBasis<4,real_t>::make(basesVec));
                }
            }
            sliceBasis->setTopology(assembler.getSlice(0)->topology());
            HeatSettings.m_STAssemblers[spaceLevels-step-1] = memory::make_shared(new gsSpaceTimesliceAssembler<real_t>(*dynamic_cast<gsSpaceTimePoissonPde<real_t>*>(SlicePde.get()),*sliceBasis,assembler.getTheta(),dirichlet::elimination,iFace::glue,true));
        }
        testMG.setUpSpaceTimeHierarchy(opt.getInt("TMG.CoarseningStrategy"));
        //   MPI_Barrier(time_comm);
        //   gsInfo<<"Rank "<<rank <<" finished ST Hierarchy.\n"<<std::flush;
        std::vector<int> spaceL = testMG.getSpaceLevels();
        std::vector<std::vector<int> > timeL(spaceLevels); for(unsigned int i=0; i<timeL.size(); i++) timeL[i].resize(0);
        for(unsigned int tLevel=0; tLevel<spaceL.size(); tLevel++) timeL[spaceL[tLevel]].push_back(tLevel);

        for(int p=0; p<gsMpi::worldSize();++p)
        {
        if(rank==p) for(unsigned int i=0; i<timeL.size(); i++)
        {
            std::cout << "rank "<<rank<<" space-level: " << i << ": ";
            for(unsigned int k=0; k<timeL[i].size(); k++)
                std::cout << timeL[i][k] << " ";
            std::cout << std::endl;
        }
        MPI_Barrier(all_comm);

        }

        gsMatrix<real_t> timeR1,timeR2;
        gsSparseMatrix<real_t> StimeR1,StimeR2;
        assembler.computeTransferForTwoSlices(0,timeR1,timeR2);
        StimeR1 = timeR1.sparseView();
        StimeR2 = timeR2.sparseView();
        //   MPI_Barrier(time_comm);
        //   gsInfo<<"Rank "<<rank <<" finished Level Setup.\n"<<std::flush;

        std::vector<real_t> h,tau;
        h.resize(timeLevels); tau.resize(timeLevels);
        std::vector<gsHeatSTSlapOperator<real_t>::Ptr> STOp(spaceLevels);


        std::vector<gsHeatSTSlapOperator<real_t>::STData > levelData(timeLevels);
        std::vector< std::vector<int> > myTimeSlices(timeLevels);
        myTimeSlices[timeLevels-1]= testMG.getMyTimeStepsIdx(timeLevels-1);
        static_cast<TPData&>(levelData[timeLevels-1]) = assembler.getTensorProductData(myTimeSlices[timeLevels-1]);
        h[timeLevels-1] = assembler.getH(0); tau[timeLevels-1] = assembler.getTau(0);

        for(int tl = timeLevels-2; tl >= 0;--tl)
        {
            myTimeSlices[tl] = testMG.getMyTimeStepsIdx(tl);
            size_t size = myTimeSlices[tl].size();
            if(spaceL[tl]==spaceL[tl+1]) //(size_t)levelData[tl+1].Kx.size()/2
            {
                levelData[tl].Kx.resize(size);
                levelData[tl].Mx.resize(size);
                levelData[tl].MWx.resize(size);
                if(useIETI)
                {
                    levelData[tl].PatchKx.resize(size);
                    levelData[tl].PatchMx.resize(size);
                    levelData[tl].PatchMWx.resize(size);
                }
                for(size_t slice = 0; slice < size;++slice)
                {
                    // No space coarsening
                    levelData[tl].Kx[slice] = levelData[tl+1].Kx[2*slice];
                    levelData[tl].Mx[slice] = levelData[tl+1].Mx[2*slice];
                    if(myTimeSlices[tl][slice]>0)levelData[tl].MWx[slice] = levelData[tl+1].MWx[2*slice];

                    if(useIETI)
                    {
                        levelData[tl].PatchKx[slice].resize(nPatches);
                        levelData[tl].PatchMx[slice].resize(nPatches);
                        if(myTimeSlices[tl][slice]>0) levelData[tl].PatchMWx[slice].resize(nPatches);
                        for(int np =0;np<nPatches;++np)
                        {
                            levelData[tl].PatchKx[slice][np] = levelData[tl+1].PatchKx[2*slice][np];
                            levelData[tl].PatchMx[slice][np] = levelData[tl+1].PatchMx[2*slice][np];
                            if(myTimeSlices[tl][slice]>0)levelData[tl].PatchMWx[slice][np] = levelData[tl+1].PatchMx[2*slice][np];
                        }

                    }

                }
                h[tl] = h[tl+1];

            }
            else
            {
                levelData[tl].Kx.resize(size);
                levelData[tl].Mx.resize(size);
                levelData[tl].MWx.resize(size);
                if(useIETI)
                {
                    levelData[tl].PatchKx.resize(size);
                    levelData[tl].PatchMx.resize(size);
                   levelData[tl].PatchMWx.resize(size);
                }
                for(size_t slice = 0; slice < size;++slice)
                {
                    // Space coarsening
                    levelData[tl].Kx[slice] = gsSparseMatrix<real_t>((prolongPtr[spaceL[tl]]->transpose())*(*levelData[tl+1].Kx[2*slice])*(*prolongPtr[spaceL[tl]])).moveToPtr();
                    levelData[tl].Mx[slice] = gsSparseMatrix<real_t>((prolongPtr[spaceL[tl]]->transpose())*(*levelData[tl+1].Mx[2*slice])*(*prolongPtr[spaceL[tl]])).moveToPtr();
                    if(myTimeSlices[tl][slice]>0)levelData[tl].MWx[slice] = gsSparseMatrix<real_t>((prolongPtr[spaceL[tl]]->transpose())*(*levelData[tl+1].MWx[2*slice])*(*prolongPtr[spaceL[tl]])).moveToPtr();

                    if(useIETI)
                    {
                        levelData[tl].PatchKx[slice].resize(nPatches);
                        levelData[tl].PatchMx[slice].resize(nPatches);
                        if(myTimeSlices[tl][slice]>0)levelData[tl].PatchMWx[slice].resize(nPatches);
                        for(int np =0;np<nPatches;++np)
                        {
                            levelData[tl].PatchKx[slice][np] = gsSparseMatrix<real_t>((PatchRestr[np][spaceL[tl]].transpose())*(*levelData[tl+1].PatchKx[2*slice][np])*(PatchRestr[np][spaceL[tl]])).moveToPtr();
                            levelData[tl].PatchMx[slice][np] = gsSparseMatrix<real_t>((PatchRestr[np][spaceL[tl]].transpose())*(*levelData[tl+1].PatchMx[2*slice][np])*(PatchRestr[np][spaceL[tl]])).moveToPtr();
                            if(myTimeSlices[tl][slice]>0)levelData[tl].PatchMWx[slice][np] = gsSparseMatrix<real_t>((PatchRestr[np][spaceL[tl]].transpose())*(*levelData[tl+1].PatchMWx[2*slice][np])*(PatchRestr[np][spaceL[tl]])).moveToPtr();
                        }

                    }
                }
                h[tl] = 2*h[tl+1];
            }
            levelData[tl].Kt.resize((size_t)levelData[tl+1].Kt.size()/2);
            levelData[tl].Mt.resize((size_t)levelData[tl+1].Kt.size()/2);
            //levelData[tl].K1t.resize((size_t)levelData[tl+1].Kt.size()/2);
            //  levelData[tl].M1t.resize((size_t)levelData[tl+1].Kt.size()/2);
            //levelData[tl].Kht.resize((size_t)levelData[tl+1].Kt.size()/2);
            //levelData[tl].Mht.resize((size_t)levelData[tl+1].Kt.size()/2);
            levelData[tl].Nt.resize((size_t)levelData[tl+1].Kt.size()/2);
            for(size_t slice = 0; slice < (size_t)(levelData[tl+1].Kt.size()/2);++slice)
            {
                int sl2 = 2*slice;// myTimeSlices[tl+1][2*slice];
                //int sl21 = 2*slice+1; // myTimeSlices[tl+1][2*slice]+1;
                int sl =  slice; // myTimeSlices[tl][slice];
                //  gsInfo<<"rank "<< rank<<" parameters: sl "<<sl<<" , sl2 "<<sl2<<" , sl21 "<<sl21<<"  total size on "<<tl<<" is: "<<(size_t)levelData[tl+1].Kt.size()/2<<" and on "<<tl+1<<" : "<<levelData[tl+1].Kt.size()<<"\n"<<std::flush;
                /*
                levelData[tl].K1t[sl] = (gsSparseMatrix<real_t>(StimeR1**levelData[tl+1].K1t[sl2]*StimeR1.transpose() + StimeR2**levelData[tl+1].K1t[sl21]*StimeR2.transpose())).moveToPtr();//levelData[tl+1].K1t[2*slice];
                levelData[tl].M1t[sl] = (gsSparseMatrix<real_t>(0.5*(StimeR1**levelData[tl+1].M1t[sl2]*StimeR1.transpose() + StimeR2**levelData[tl+1].M1t[sl21]*StimeR2.transpose()))).moveToPtr();//levelData[tl+1].M1t[2*slice];
                levelData[tl].Kht[sl] = (gsSparseMatrix<real_t>((StimeR1**levelData[tl+1].Kht[sl2]*StimeR1.transpose() + StimeR2**levelData[tl+1].Kht[sl21]*StimeR2.transpose()))).moveToPtr();//levelData[tl+1].Kht[2*slice];
                levelData[tl].Mht[sl] = (gsSparseMatrix<real_t>(0.5*(StimeR1**levelData[tl+1].Mht[sl2]*StimeR1.transpose() + StimeR2**levelData[tl+1].Mht[sl21]*StimeR2.transpose()))).moveToPtr();//levelData[tl+1].Mht[2*slice];
                levelData[tl].Kt[sl] = (gsSparseMatrix<real_t>(*levelData[tl].K1t[sl] + (*levelData[tl].Kht[sl]))).moveToPtr() ;
                levelData[tl].Mt[sl] = (gsSparseMatrix<real_t>(*levelData[tl].M1t[sl] + (*levelData[tl].Mht[sl]))).moveToPtr() ;
*/
                //if(myTimeSlices[tl][slice]>0)levelData[tl].Nt[sl] = levelData[tl+1].Nt[sl2];

                //   if(sl==0)gsInfo<<"Kt: \n"<<levelData[tl+1].Kt[sl2]->toDense()<<"\n\n";
                //   if(sl==0)gsInfo<<"Mt: \n"<<levelData[tl+1].Mt[sl2]->toDense()<<"\n\n";
                //   if(sl==0)gsInfo<<"timeR1: \n"<<timeR1<<"\n\n";
                //   if(sl==0)gsInfo<<"timeR2: \n"<<timeR2<<"\n\n";

                //  levelData[tl].M1t[sl] = (gsSparseMatrix<real_t>(2*(*levelData[tl+1].M1t[sl2]))).moveToPtr();

                levelData[tl].Mt[sl] = (gsSparseMatrix<real_t>(2*(*levelData[tl+1].Mt[sl2]))).moveToPtr();
                levelData[tl].Kt[sl] = levelData[tl+1].Kt[sl2];

                if(sl>0)levelData[tl].Nt[sl] = levelData[tl+1].Nt[sl2];
            }
            tau[tl]  =tau[tl+1]*2;

        }

        for(int p=0; p<gsMpi::worldSize();++p)
        {
            if(detailedPrint && p==rank)
            {
                gsInfo<<"TimeSlices for rank "<<rank<<": \n";
                for(size_t i=0; i< myTimeSlices.size();++i)
                {
                    gsInfo<<"Level: "<<i<<" ";
                    for(size_t j=0; j<myTimeSlices[i].size();++j) gsInfo<<myTimeSlices[i][j]<<", ";
                    gsInfo<<"\n";
                }
                gsInfo<<"\n";
            }
            MPI_Barrier(all_comm);
        }

        if(rank==0) for(int tl = 0; tl<timeLevels;++tl)
        {
            gsInfo << "Level: "<<tl<<" - h: "<<h[tl]<<" -tau: "<<tau[tl]<<" - tau/h^2: "<< tau[tl]/(h[tl]*h[tl])<<"\n"<<std::flush;
        }

        //      MPI_Barrier(time_comm);
        //   gsInfo<<"Rank "<<rank <<" Building Matrix hierachy finished.\n"<<std::flush;
        if(rank==0)gsInfo<<"Building Matrix hierachy finished!\n"<<std::flush;
        time.restart();
        for(int step=0; step<spaceLevels; step++)
        {
            std::vector<gsHeatSTSlapOperator<real_t>::STData> data(timeL[step].size());
            for(size_t i=0; i<timeL[step].size();++i) data[i] = levelData[timeL[step][i]];
            STOp[step] =  gsHeatSTSlapOperator<real_t>::make(timeL[step], timeLevels,spaceLevels,step,myTimeSlices, data,prolongPtr, timeR1,timeR2, tau[timeLevels-1], math::exp2(spaceLevels-1-step)*h[timeLevels-1], HeatSettings, rank==1);
            STOp[step]->setOutputFilename(opt.getString("iFilename"));
            STOp[step]->init();
            testMG.setSTSlapOperator(STOp[step], step);
        }
        real_t initTime = time.stop();
        if(rank==0)  std::cout <<  std::endl << "Init MG -- finished in: " << initTime << " sec." <<  std::endl<<std::flush;

        //   MPI_Barrier(time_comm);
        //    gsInfo<<"Rank "<<rank <<" Operator hierachy finished.\n"<<std::flush;
        if(rank==0) testMG.printSpaceTimeHierarchy();

        if(rank==0) std::cout << "init sol vectors" << std::endl;
        MPI_Barrier(all_comm);

        STMultigridMPI<gsHeatSTSlapOperator<real_t> >::STSolVector vecT;
        testMG.initSTSolVectors(vecT, false);
        std::vector<gsMatrix<real_t> > flT = levelData[timeLevels-1].rhs;

        if(rank==0)  std::cout << "start ST-MG..." <<  std::endl<<std::flush;
        time.restart();

        int nIt=0;
        switch(opt.getInt("SolverType"))
        {
        case 0:
            testMG.STForwardSolve(vecT,flT,timeLevels-1);
            break;
        case 1:
            nIt=testMG.TimeMGM_MPI(vecT, flT, opt.getInt("TMG.MaxIters"), opt.getReal("TMG.Tolerance"),-1, (rank==0) | detailedPrint, true);
            break;
        case 2:
            nIt=testMG.PGMRESSolve(vecT,flT, opt.getInt("TMG.MaxIters"), opt.getReal("TMG.Tolerance"),-1, (rank==0) | detailedPrint, true);
            break;
        }

        real_t solveTime = time.stop();
        if(rank==0)  std::cout <<  std::endl << "finished in: " << solveTime << " sec." <<  std::endl<<std::flush;

        double res = testMG.ComputeSpaceTimeResiduum(vecT, flT, timeLevels-1);
        if(rank==0)  std::cout << "\tResiddum: " << res <<  std::endl<<std::flush;
        if(rank==0)     gsInfo<<"Dofs: "<<assembler.numDofs()<<"\n";

        MPI_Barrier(all_comm);
        solVMG.setZero(assembler.numDofs(),1);
        if(rank==0)     gsInfo<<"start building global solution\n";
        testMG.buildGlobalSolution(vecT.u[timeLevels-1], assembler.getBlockViewSizes(), solVMG);
        if(rank==0)     gsInfo<<"finished global solution\n";
        if(rank==0)
        {

            // gsVector<index_t> sizes = assembler.getBlockViewSizes();
            //  gsVector<index_t> cols(1);
            //  cols<<1;
            // gsMatrix<real_t>::BlockView view = solVMG.blockView(sizes,cols);
            // for(size_t slice = 0; slice <view.numRowBlocks();++slice)
            //     view(slice,0) = vecT.u[timeLevels-1][slice];

            //  gsInfo<<"SolV of MG\n: "<<solVMG.transpose()<<"\n\n";
       //     Permutation perm = assembler.getTensorPermutation();
         //   solVMG = perm.inverse()*solVMG;
            gsField<> solMG = assembler.constructSolution(solVMG);


            if (opt.getSwitch("p"))
            {
                gsField<real_t> field = gsFieldCreator<real_t>::absError(solMG,exact);
                gsWriteParaview<>( field, "ST_MGdiff"+util::to_string(0), 1000);
                gsWriteParaview<>( solMG, "ST_SolMG", 1000,true);
                gsFileManager::open("ST_SolMG.pvd");
            }

        }

        if(rank==0 && opt.getString("filename").compare("")!=0)
        {
            std::fstream out( opt.getString("filename"),std::fstream::app|std::fstream::out);
            int maxXDegree = 0;
            for(int d=0; d<pde.domain().parDim()-1;++d)
                if(maxXDegree < refine_bases.maxDegree(d)) maxXDegree = refine_bases.maxDegree(d);
            int maxTDegree = refine_bases.maxDegree(pde.domain().parDim()-1);

            out<<"#dofs: "<<assembler.numDofs()<<", ref (x,t): ("<<opt.getInt("refX")<<","<<opt.getInt("refT")<<"), deg (x,t): ("<<maxXDegree<<","<<maxTDegree<<"), #slaps: "<<
                 nTimeSlices<<", #t-levels: "<<timeLevels<<",\t #it: "<<nIt<<", residuum: "<<res <<", time: Assemble-Init-Solve: "<<assTime<<", "<<initTime<<", "<<solveTime<<"\n";

        }

    }


    if(false && rank==0)
    {
        // For comparism, the exact solution
        assembler.assemble();
        assembler.printTimings();
     //   Permutation perm = assembler.getTensorPermutation();
        gsSparseMatrix<> K = assembler.matrix();
    //    K=K.twistedBy(perm);


        gsMatrix<real_t> f = assembler.rhs();
     // gsMatrix<real_t> f = perm*assembler.rhs();
        //  gsMatrix<real_t> mat1,mat2;
        //   assembler.computeTransferForTwoSlices(0,mat1,mat2);

        // gsInfo<<"K: "<<K.toDense()<<"\n";
        //  gsInfo<<"f: "<<f.transpose()<<"\n";

#if defined(GISMO_WITH_PARDISO)
        gsSparseSolver<>::PardisoLU solver;
        //gsSparseSolver<>::PardisoLU::ParameterType& param =  solver.pardisoParameterArray();

        // param[10] = 0;
        //  param[12] = 0;
#elif defined(GISMO_WITH_SUPERLU)
        gsSparseSolver<>::SuperLU solver;
#else
        gsSparseSolver<>::LU solver;
#endif

        solver.compute(K);
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eiK(K.toDense());
        //gsInfo<<"eigenvalues K: "<<eiK.eigenvalues().transpose()<<"\n";
        gsMatrix<> solV = solver.solve( f );
        // gsInfo<<"solution:\n"<<solV.transpose()<<"\n";
   //     solV = perm.inverse()*solV;
        gsField<> sol = assembler.constructSolution(solV);
        gsInfo<<"Dofs: "<<assembler.numDofs()<<"\n";
        //   gsFunctionExpr<> gder("pi*cos(pi*(x))*sin(pi*(y))","pi*sin(pi*x)*cos(pi*y)",2);
        gsSpaceTimeSliceNorm<real_t> norm(sol,exact,assembler.getdGInterfaces(), assembler.getInitialBoundary(),assembler.getTopBoundary());
        norm.setTheta(assembler.getTheta());
        real_t l2error = sol.distanceL2(exact);
        gsInfo<<"l2 error ||uhMG - uh||: "<<(solV - solVMG).norm()<<"\n";
        gsInfo<<"L2 error ||uh - g||: "<<l2error<<"\n";
        gsInfo<<"dG error ||uh - g||: "<<norm.compute()<<"\n";

        // Plot solution in paraview

        if (opt.getSwitch("p"))
        {
            gsWriteParaview<>( sol, "ST_Sol", 1000,true);

            // Write approximate and exact solution to paraview files
            gsInfo<<"Plotting in Paraview...\n";
            const gsField<> exactF( pde.domain(), exact, false );
            const gsField<> rhs( pde.domain(), *pde.rhs(), false );
            gsWriteParaview<>( exactF, "ST_Ex", 1000);
            gsWriteParaview<>( rhs, "ST_rhs", 1000);

            // Run paraview
            gsFileManager::open("ST_Sol.pvd");


            gsField<real_t> field = gsFieldCreator<real_t>::absError(sol,exact);
            gsWriteParaview<>( field, "ST_diff"+util::to_string(0), 1000);

            // Run paraview
            //     gsFileManager::open(("paraview ST_diff"+util::to_string(0)+".pvd  &").c_str());

        }
    }
}



//Important Note: The boundary condition has to be called with
// bcInfo.addCondition(***, *** , *** , g);
// in order that the parallelization works correctly! So no &g !!!
// Otherwise, you have a data race!


int main (int argc, char** args)
{
    gsMpi& mpi = gsMpi::init(argc, args);
    gsMpiComm comm(mpi.worldComm());

    std::string filename="";
    std::string iFilename="";

    int testcase =1;
    bool plot = false;
    int test=0;
    int spaceLevel = 1;
    int nSlices = 2;
    real_t Tlen = 1;
    real_t theta = 1;

    int numRefineX  = 2;
    int numRefineT  = 1;
    int numElevateX = 0;
    int numElevateT = 0;

    real_t TMGtol = 1.e-8;
    int TMGmaxIter = 100;
    int Presmoothing = 2;
    int Postsmoothing = 2;
    int cycles = 1;
    real_t omega = 0.5;
    bool hybrid = false;
    bool tridiag = false;
    bool detailedOutput = false;
    int cStrategy = 4;
    int sStrategy = 1;

    int solver = 0;
    int maxIter = 10;
    real_t tol = 1.e-3;
    real_t PTol = 1.e-4;
    int PIt = 40;
    int decompostion = 2;
    int decompostionSolv = 0;
    int NCoresTime = 1;
    int NCoresSpace = 1;
    int n = 2;

    std::string m = "C";


    gsCmdLine cmd("Hi, I will test IETIdG");
    cmd.addSwitch("p","plot", "Plot result in ParaView format",plot);
    cmd.addString("", "iFilename","filename for data output", iFilename);
    cmd.addString("", "filename","filename for data output", filename);
    cmd.addInt("c","case","Choosen testcase",testcase);
    cmd.addInt("n", "nPatches", "number of spatial patches in one dimension",n);
    cmd.addInt("","test","(0) experiment, (1) convergence test", test);

    cmd.addInt("","refX","Number of refinements",numRefineX);
    cmd.addInt("","refT","Number of refinements",numRefineT);
    cmd.addInt("","elX","maximal elevation in spatial direction",numElevateX);
    cmd.addInt("","elT","maximal elevation in temporal direction",numElevateT);
    cmd.addReal("", "theta", "theta parameter for the dG-scheme",theta);
    cmd.addReal("", "Tlen", "length of a Timeslab",Tlen);
    cmd.addSwitch("","detailedOutput","prints more stuff, may be huge for large number of procs!", detailedOutput);
    cmd.addInt("","SolverType","Solver Type for the ST-system: (0) slapwise-direct, (1) MG-iteration, (2) MG-GMRES", sStrategy);

    cmd.addInt("t","nSlices", "number of time slices", nSlices );
    cmd.addInt("", "MG.Levels", "number of spaceLevels for MG", spaceLevel);

    cmd.addReal("", "TMG.Tolerance", "tolerance for the space time MG", TMGtol);
    cmd.addInt("", "TMG.MaxIters", "maximal number of Iteration", TMGmaxIter);
    cmd.addInt("", "TMG.NumPreSmooth", "number of smoothing steps", Presmoothing);
    cmd.addInt("", "TMG.NumPostSmooth", "number of smoothing steps", Postsmoothing);
    cmd.addReal("", "TMG.Damping", "damping parameter omega", omega);
    cmd.addInt("", "TMG.NumCycles", "number of cycles", cycles);
    cmd.addSwitch("", "TMG.TriDiag", "handle TriDiag formulation", tridiag);
    cmd.addSwitch("","TMG.Hybrid","use hybrid smoother", hybrid);
    cmd.addInt("","TMG.CoarseningStrategy","the coarsening strategy: (0) no space, (1) low space, (2) moderate space, (3) aggressive space, (4) automatic (default), (5) stochastic",cStrategy);
    cmd.addInt("", "TimeCores","Number of cores used for the time paralellization", NCoresTime);
    cmd.addInt("", "SpaceCores", "Cores for space parallelization.",NCoresSpace);

    cmd.addInt("", "TMGSlap.Solver", "(0) direct, (1) Slicewise Direct", solver);
    cmd.addInt("", "TMGSlap.Iterations", "maximum number of iterations", maxIter);
    cmd.addInt("", "TMGSlap.DecompositionMethod","Decomposition method for time direction: (0) diagonalization (may not work), (1) Real Schur, (2) Complex Schur", decompostion);
    cmd.addInt("", "TMGSlap.DecompositionSolver","Solver for the remaining spatial problems: (0) direct, (1) IETI", decompostionSolv);
    cmd.addReal("", "TMGSlap.Tolerance", "tolerance for the iterative solver for each subproblem of A", tol );
    cmd.addInt("","TMGSlap.PreconditionIterations", "max number of iteraions for preconditioner in approximative solve",PIt);
    cmd.addReal("","TMGSlap.PreconditionTol", "tolerance for the preconditioner in approximative solve", PTol);
    cmd.addString("s", "TMGSlap.IETI.Strategy", "Choosen strategy",m);

    try { cmd.getValues(argc,args); } catch (int rv) { return rv; }

    gsOptionList option = cmd.getOptionList();
    if(mpi.worldRank()==0) option.print(gsInfo);

    //nitsche = true;


    int dim;
    if(testcase>0 && testcase<5) dim = 2;
    if(testcase>=5 && testcase <9) dim = 3;
    if(testcase>=9) dim = 4;

    if(dim == 4)
        plot  = false;
    if(mpi.worldRank()==0) gsInfo<<"Dimension of Domain: "<<dim<<"\n";

    //////// Right-hand side and analytical solution ////////
    // Define source function
    gsFunctionExpr<> f,f_x,f_t,g;

    if(dim==2)
    {
        f= gsFunctionExpr<>("pi*sin(pi*x)*(cos(pi*y)+pi*sin(pi*y)) ",dim);
        f_x= gsFunctionExpr<>("pi*sin(pi*x)",dim-1);
        f_t=gsFunctionExpr<>("cos(pi*x)+pi*sin(pi*x)",1);
        g=gsFunctionExpr<>("sin(pi*x)*sin(pi*y)+1",dim);

    }
    else if(dim==3)
    {
        f= gsFunctionExpr<>("pi*sin(pi*x+1)*sin(pi*y)*(cos(pi*z+2)+2*pi*sin(pi*z+2)) ",dim);
        f_x= gsFunctionExpr<>("pi*sin(pi*x+1)*sin(pi*y)",dim-1);
        f_t=gsFunctionExpr<>("cos(pi*x+2)+2*pi*sin(pi*x+2)",1);
        g=gsFunctionExpr<>("sin(pi*x+1)*sin(pi*y)*sin(pi*z+2)",dim);

        //  f= gsFunctionExpr<>("if((x-0.7)^2+(y-4.5)^2<0.04,5,0)",dim);
        //   f_x= gsFunctionExpr<>("if((x-0.4)^2+(y-4.6)^2<0.04,5,0)",dim-1);
        //    f= gsFunctionExpr<>("0",dim);
        //    f_x= gsFunctionExpr<>("0",dim-1);
        //   f_t=gsFunctionExpr<>("0",1);
        //   g=gsFunctionExpr<>("10",dim);
    }
    else if(dim==4)
    {
        f= gsFunctionExpr<>("pi*sin(pi*x)*sin(pi*y)*sin(pi*z)*(cos(pi*w)+3*pi*sin(pi*w)) ",dim);
        f_x= gsFunctionExpr<>("pi*sin(pi*x)*sin(pi*y)*sin(pi*z)",dim-1);
        f_t=gsFunctionExpr<>("cos(pi*x)+3*pi*sin(pi*x)",1);
        g=gsFunctionExpr<>("sin(pi*x)*sin(pi*y)*sin(pi*z)*sin(pi*w)",dim);
    }

    //  gsFunctionExpr<> dn_g("pi*sin(pi*x)",dim-1);//Neumann at x=1
    // Print out source function and solution
    if(mpi.worldRank()==0) gsInfo<<"Source function "<< f << "\n";
    if(mpi.worldRank()==0) gsInfo<<"Exact solution "<< g <<"\n" << "\n";

    /////////////////// Setup geometry ///////////////////
    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches;


    std::vector<patchSide> connectedComponent;
    switch(testcase){
    case 1:
        patches = gsNurbsCreator<>::BSplineSquareGrid(1,nSlices, 1);
        break;
    case 2:
        patches = gsNurbsCreator<>::BSplineSquareGrid(2,nSlices, 0.5);
        break;
    case 3:
        patches = gsNurbsCreator<>::BSplineSquareGrid(3,nSlices, 1);
        break;
    case 4:
        patches = gsNurbsCreator<>::BSplineSquareGrid(4,nSlices, 0.5);
        break;
    case 5:
        patches = gsNurbsCreator<>::BSplineCubeGrid(1,1,nSlices,1);
        break;
    case 6:
        patches = gsNurbsCreator<>::BSplineCubeGrid(2,1,nSlices,0.5,+0.25,+0.25);
        break;
    case 7:
        patches = gsNurbsCreator<>::BSplineCubeGrid(n,n,nSlices,1);
        break;
    case 8:
    {
	/*
        gsFileData<> fileData("yeti3d_u21p.xml");
        if (fileData.getFirst< gsMultiPatch<> >(patches))
        {
            int nPatches= patches.nPatches();
            for(int slice = 1; slice<nSlices;++slice)
            {
                for(int np = (slice-1)*nPatches; np< slice*nPatches;++np)
                {
                    gsGeometry<real_t>::uPtr newPatch = patches.patch(np).clone();
                    gsMatrix<real_t>& coefs = newPatch->coefs();
                    real_t diff = coefs.col(2).maxCoeff() - coefs.col(2).minCoeff();
                    for(int j=0; j<coefs.rows();++j)
                        coefs.col(2)(j)+=diff;
                    patches.addPatch( give(newPatch) );
                }
            }
            for(index_t np =0; np<patches.nPatches();++np)
                patches.patch(np).degreeElevate(1,2);
            patches.computeTopology();
        }

        else
            return 1;
        break;
        //*/

        gsFileData<> fileData_temp ( "yeti_mp2.xml"); //2d version
        gsMultiPatch<> ieti;
        if (fileData_temp.getFirst<gsMultiPatch<> >(ieti))
        {

           // patches = patches.uniformSplit();
           // ieti = ieti.uniformSplit();
            ieti.computeTopology();
            patches.clear();
            for(size_t i =0; i<ieti.nPatches();i++)
            {
                gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&ieti.patch(i));
                patches.addPatch(gsNurbsCreator<>::lift3D(*tb));
            }

            int nPatches= patches.nPatches();
            connectedComponent = getConnectedBoundaryComponent(ieti,*ieti.bBegin()); //find the outer boundary
            for(int slice = 1; slice<nSlices;++slice)
            {
                for(int np = (slice-1)*nPatches; np< slice*nPatches;++np)
                {
                    gsGeometry<real_t>::uPtr newPatch = patches.patch(np).clone();
                    gsMatrix<real_t>& coefs = newPatch->coefs();
                    real_t diff = coefs.col(2).maxCoeff() - coefs.col(2).minCoeff();
                    for(int j=0; j<coefs.rows();++j)
                        coefs.col(2)(j)+=diff;
                    patches.addPatch(give(newPatch));
                }
            }
            for(size_t np =0; np<patches.nPatches();++np)
                patches.patch(np).degreeElevate(1,2);
            patches.computeTopology();
        }

        else
            return 1;
        break;
        //*/
    }
    case 10:
    {
        patches = gsNurbsCreator<>::BSplineCubeGrid(2,2,2,1);
        gsMultiPatch<real_t> patches4D;
        for(size_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<3,real_t>* tb = static_cast<gsTensorBSpline<3,real_t>* >(&patches.patch(i));
            patches4D.addPatch(gsNurbsCreator<>::lift4D(*tb));
        }

        int nPatches= patches4D.nPatches();
        for(int slice = 1; slice<nSlices;++slice)
        {
            for(int np = (slice-1)*nPatches; np< slice*nPatches;++np)
            {
                gsGeometry<real_t>::uPtr newPatch = patches4D.patch(np).clone();
                gsMatrix<real_t>& coefs = newPatch->coefs();
                real_t diff = coefs.col(3).maxCoeff() - coefs.col(3).minCoeff();
                for(int j=0; j<coefs.rows();++j)
                    coefs.col(3)(j)+=diff;
                patches4D.addPatch(give(newPatch));
            }
        }
        patches4D.computeTopology();
        patches = give(patches4D);

        break;
    }
    default:
        if(mpi.worldRank()==0) gsInfo<<"No valid geometry given! \n";
        return 1;
        break;
    }

    //scale the length of the timeslap to len
    real_t len = patches.patch(0).coefs().col(dim-1).maxCoeff() - patches.patch(0).coefs().col(dim-1).minCoeff();
    for(size_t np = 0; np<patches.nPatches();++np)
        patches.patch(np).scale(Tlen/len,dim-1);

    /*
#ifdef _OPENMP
    // To get the number of threads from the environment
    int nThreads;
    const char * nProcs = getenv("OMP_NUM_THREADS");
    if(nProcs != NULL)
        sscanf( nProcs, "%d", &nThreads );
    else
        nThreads = 1; //one OpenMP thread
    nThreads = math::min(nThreads,(int)patches->nPatches());

    std::ostringstream oss;
    oss << nThreads;

    setenv("OMP_NUM_THREADS", oss.str().c_str(),1);
#endif
*/

    /////////////////// Setup boundary conditions ///////////////////
    // Define Boundary conditions
    gsBoundaryConditions<> bcInfo;

    switch(testcase){
    case 1:
    case 2:
    case 3:
    case 4:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
        {
            if(it->side() == boundary::south)
            {
                bcInfo.add(it->patch, boundary::south, "Initial", g);
                // bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, f);
            }
            else if(it->side() != boundary::north)
                bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, g);
        }
        break;
    case 5:
    case 6:
    case 7:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
        {
            if(it->side() == boundary::front)
            {
                bcInfo.add(it->patch, boundary::front, "Initial", g);
                // bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, f);
            }
            else if(it->side() != boundary::back && it->side() == boundary::right)
                bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, g);
        }
        break;
    case 8:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
        {
            if(it->side() == boundary::front)
            {
                bcInfo.add(it->patch, boundary::front, "Initial", g);
                // bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, f);
            }
          //  else if(it->side() != boundary::back)
           //     bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, g);
        }

        for(int slice = 0; slice<nSlices;++slice)
        {
            for(std::vector<patchSide>::iterator it=connectedComponent.begin();it!=connectedComponent.end();++it)
            {
                bcInfo.addCondition(slice*(patches.nPatches()/nSlices)+it->patch, it->side(), condition_type::dirichlet, g);
            }

        }
        //*/
        break;
    case 9:
        for(gsMultiPatch<real_t>::const_biterator it=patches.bBegin();it!=patches.bEnd();it++)
        {
            if(it->side() == boundary::stime)
            {
                bcInfo.add(it->patch, it->side(), "Initial", g);
                //   bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, g);
            }
            else if(it->side() != boundary::etime)
                bcInfo.addCondition(it->patch,it->side(), condition_type::dirichlet, g);
        }
        break;
    default:
        if(mpi.worldRank()==0) gsInfo<<"No valid geometry given! \n";
        return 1;
        break;
    }
    /*
    if(test3D)
    {
        gsMultiPatch<real_t> patches3D;
        for(index_t i =0; i<patches.nPatches();i++)
        {
            gsTensorBSpline<2,real_t>* tb = static_cast<gsTensorBSpline<2,real_t>* >(&patches.patch(i));
            patches3D.addPatch(gsNurbsCreator<>::lift3D(*tb));
            bcInfo.add(i, boundary::front, "Initial", g);
            if(testcase == 6)
                patches3D.patch(i).degreeElevate(1,2);

        }
        patches3D.computeTopology();

        patches = give(patches3D);

    }
    */
    if (mpi.worldRank()==0 && plot)
    {
        gsWriteParaview<>( patches, "IETI_patch", 1000 );
        gsFileManager::open("IETI_patch.pvd");
    }

    ////////////////////// Refinement h and p //////////////////////
    // Refinement

    // Copy basis from the geometry
    gsMultiBasis<> refine_bases( patches );

    // Find maximum degree with respect to all the variables
    refine_bases.setDegree(refine_bases.maxCwiseDegree());
    for(int d=0; d<dim-1;++d)
        refine_bases.degreeIncrease(numElevateX,d);
    refine_bases.degreeIncrease(numElevateT,dim-1);

    int maxXDegree = 0;
    for(int d=0; d<dim-1;++d)
        if(maxXDegree < refine_bases.maxDegree(d)) maxXDegree = refine_bases.maxDegree(d);
    int maxTDegree = refine_bases.maxDegree(dim-1);

    if(mpi.worldRank()==0) gsInfo<<"Degree in X set to: "<<maxXDegree<<"\n";
    if(mpi.worldRank()==0) gsInfo<<"Degree in T set to: "<<maxTDegree<<"\n";

    if(false && (testcase == 5 || testcase == 7)) //Not supported for MG.
    {
        for(size_t i = 0; i < patches.nPatches(); ++i)
        {
            if(i%2 == 0)
                refine_bases[i].uniformRefine();
        }
    }
    std::fstream out(iFilename,std::ofstream::out | std::ofstream::app);

    out<<"deg (x,t): "<<maxXDegree<<","<<maxTDegree<<"\t ref(x,t): "<<option.getInt("refX")<<","<<option.getInt("refT")<<"\n"<<std::flush;
    ///////////////////////ASSEMBLE////////////////////////////////////

    if(mpi.worldRank()==0) gsInfo<<bcInfo<<"\n";


    gsPiecewiseFunction<real_t>  alpha;
    gsFunctionExpr<> a1("1.e+0",dim-1);
    gsFunctionExpr<> a2("1.e+0",dim-1);
    for(size_t np = 0; np<patches.nPatches();np++)
        if(np%2 == 0)
            alpha.addPiece(a1);
        else
            alpha.addPiece(a2);


    gsSpaceTimePoissonPde<real_t> pde(patches,bcInfo,f,f_x,f_t,alpha,&g);

    switch(test)
    {
    case 0:
        doExperiment(pde, refine_bases,option,g);
        break;
    }





    //return  1 for failures and 0 for success
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

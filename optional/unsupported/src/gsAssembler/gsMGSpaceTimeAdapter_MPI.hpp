/**  gsMGSpaceTimeAdapter_MPI.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    Created on:  2017-11-23
*/


#include  <gsAssembler/gsMGSpaceTimeAdapter_MPI.h>
#include <gsIETI/gsIETIAdapter.h>
#include <gsSolver/gsSumOp.h>
#include <gsIETI/gsDistributedBlockOp.h>
#include <gsIETI/gsDistributedScaledOp.h>
#include <gsIETI/gsDistributedProductOp.h>
#include <gsIETI/gsDistributedSumOp.h>
#include <gsIETI/gsParallelMinRes.h>
#include <gsIETI/gsParallelCG.h>
#include <gsSolver/gsParallelPreconditioner.h>
#include <fstream>
#include <gsSolver/gsComplexify.h>

#include <gsCore/gsConnectedBoundaryComponent.h>
#include <gsMultiGrid/gsParallelGridHierarchy.h>
#include <gsMultiGrid/gsParallelMultiPatchPreconditioners.h>
#include <gsMultiGrid/gsMassSmoother.h>
#include <gsAssembler/gsParameterDomainAssembler.h>
#include <gsSolver/gsIterativeSolverOp.h>

#ifdef GISMO_WITH_MPI
namespace gismo {


template<typename T>
gsSpaceParallelHSTSlapOp<T>::gsSpaceParallelHSTSlapOp(real_t h, real_t tau, int NTimeLevels) :Base(h,tau,NTimeLevels) {}

//Matrices are already in Tensor form
template<typename T>
gsSpaceParallelHSTSlapOp<T>::gsSpaceParallelHSTSlapOp(const std::vector<int> &TimeLevels, int NTimeLevels,
                                                      int spaceLevels,
                                                      int spaceLevel,
                                                      std::vector<std::vector<int> > myTimeSlices,
                                                      const gsSortedVector<size_t>& myPatches,
                                                      std::vector<STData > data,
                                                      std::vector<typename gsParallelOperator<T>::Ptr > prolongation,
                                                      std::vector<typename gsParallelOperator<T>::Ptr > restriction,
                                                      std::vector<typename gsPatchSubassembledTopology<T>::Ptr>  subassTopology,
                                                      std::vector<typename gsParallelGlobalLocalHandler::Ptr> handlers,
                                                      std::vector<typename gsParallelOperator<T>::Ptr> Kx,
                                                      std::vector<typename gsParallelOperator<T>::Ptr> Mx,
                                                      std::vector<typename gsParallelOperator<T>::Ptr> Mwx,
                                                      gsMatrix<T> timeRestr1,
                                                      gsMatrix<T> timeRestr2,
                                                      real_t tau,
                                                      real_t h,
                                                      typename Base::HeatSTSlapSettings settings, const gsMpiComm& comm, bool output )
    :Base(TimeLevels,NTimeLevels,spaceLevels,spaceLevel,myTimeSlices,data,timeRestr1,timeRestr2,tau,h,settings,output),
      myPatches_(myPatches),comm_(comm), m_subassTopology(subassTopology), m_parSpaceHandler(handlers), m_prolongation(prolongation),m_restriction(restriction),
      m_Kx(Kx), m_Mx(Mx), m_MWx(Mwx)
{
    if(TimeLevelsInp_.size()==0) return;

    int nRestr = restriction.size();
    if(output_)std::cout<<"finished restr "<< std::endl;
    transferSpace.resize(nRestr);
    transferTSpace.resize(nRestr);

    for(unsigned int step=0; step<transferSpace.size(); step++)
    {
        transferTSpace[step] = gsKroneckerOp<T>::make(gsIdentityOp<T>::make(timeRestr1_.cols()),prolongation[step]);
        transferSpace[step] = gsKroneckerOp<T>::make(gsIdentityOp<T>::make(timeRestr1_.cols()),restriction[step]);
    }
    if(output_)std::cout<<"finished Space transfer "<< std::endl;


}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::init()
{
    bool outputIETI_iter = true;
    real_t tolExactInnerSolv = 1.e-10;
    real_t complexBound = 1.e-10;

    if(output_)std::cout<<" Start initializing "<< std::endl;
    A.resize(TimeLevelsInp_.size());
    B.resize(TimeLevelsInp_.size());


    const gsMultiPatch<T>&  mp = settings_.m_STAssemblers[spaceLevel]->patches();
    const gsMultiBasis<T>&  bases = settings_.m_STAssemblers[spaceLevel]->multiBasis();
    const gsDofMapper & globalMapper = settings_.m_STAssemblers[spaceLevel]->system().colMapper(0);
    gsAssembler<T>* tempAss = settings_.m_STAssemblers[spaceLevel]->clone();

   // gsInfo<<"bases: "<<bases.totalSize()<<"\t"<<"mapper: "<<globalMapper.freeSize()<<"\n\n";

    std::vector<gsDofMapper> locMappers(mp.nPatches());
    for(size_t i=0; i<mp.nPatches();++i)
    {
        gsPde<T>* pde = settings_.m_STAssemblers[spaceLevel]->pde().restrictToPatch(i);
        tempAss->initialize(*pde,bases.basis(i),settings_.m_STAssemblers[spaceLevel]->options());

        locMappers[i] = tempAss->system().colMapper(0);
    }
    delete tempAss;

    //gsInfo<<"GlobalMapper:\n";
   // globalMapper.print(gsInfo);

    typename gsPatchSubassembledTopology<T>::Ptr top = gsPatchSubassembledTopology<T>::make(myPatches_,bases,locMappers,globalMapper);
    top->reorderLike(globalMapper);
    gsDofMapper procLocalPatchMapper= top->getLocalPatchMapper();

   // gsInfo<<"procLocalPatchMapper :\n";
  //  procLocalPatchMapper.print(gsInfo);

    std::vector<gsDofMapper> procLocalMapper(comm_.size());
    procLocalMapper[comm_.rank()] = top->getLocalMapper();

    gsPatchInterfaceConnections<T> patchConnection(myPatches_,bases,locMappers,procLocalPatchMapper,globalMapper,comm_);
    patchConnection.init();
    gsDofMapper procGlobMapper = patchConnection.generateProcGlobalMapper();

    m_parHandler =  gsParallelGlobalLocalHandler::make(procGlobMapper,procLocalMapper,comm_);




    for(unsigned int level=0; level<TimeLevelsInp_.size(); level++)
    {
        if(output_)std::cout<<" Do level "<<level<<"."<< std::endl<<std::flush;
        int size = (int)myTimeSlices[TimeLevelsInp_[level]].size();

        A[level].resize(size);
        B[level].resize(size);
        for(int slice = 0; slice<size;++slice) //data[level].Kt.size()
        {
            int tSlice = myTimeSlices[TimeLevelsInp_[level]][slice];
            size_t nPatches = settings_.m_spatialAssemblers[spaceLevel]->patches().nPatches();
            //   if(output_)std::cout<<" , level "<<level<<", slice "<<slice <<"  calculating A,B, mapped slice "<<tSlice <<", nPatches: "<<nPatches<< std::endl;


            std::vector<typename gsLinearOperator<T>::Ptr> localMatricesA(nPatches),localMatricesB(nPatches);

            for(size_t i=0; i< myPatches_.size();++i)
            {
                size_t np = myPatches_[i];
                //  if(output_)std::cout<<"A: Do the "<<i<<"th patch ("<<np<<")" << std::endl<<"size data: "<<data_.size()<<" \t size data.Mat: "<<data_[level].Kt.size()<<" and " <<data_[level].PatchMx.size()<<"\n";
                //  if(output_) std::cout<<"Op sizes: "<<data_[level].Kt[tSlice]->rows()<<" - "<<data_[level].PatchMx[slice][np]->rows()<<" - "<<data_[level].Mt[tSlice]->rows()<<" - "<< data_[level].PatchKx[slice][np]->rows()<<"\n";
                localMatricesA[np] = gsSumOp<T>::make(gsKroneckerOp<T>::make(makeMatrixOp(data_[level].Kt[tSlice]),makeMatrixOp(data_[level].PatchMx[slice][np])),
                                                      gsKroneckerOp<T>::make(makeMatrixOp(data_[level].Mt[tSlice]),makeMatrixOp(data_[level].PatchKx[slice][np])));

                if(tSlice>0)
                {
                    //      if(output_)std::cout<<"B Do the "<<i<<"th patch ("<<np<<")" << std::endl<<"size data: "<<data_.size()<<" \t size data.Mat: "<<data_[level].Nt.size()<<" and " <<data_[level].PatchMWx.size()<<"\n";
                    //      if(output_) std::cout<<"PatchWMx size: "<<data_[level].PatchMWx[slice].size()<<"\n";
                    //      if(output_) std::cout<<"Op sizes: "<<data_[level].Nt[tSlice]->rows()<<" - "<<data_[level].PatchMWx[slice][np]->rows()<<"\n";
                    localMatricesB[np] = gsKroneckerOp<T>::make(makeMatrixOp(data_[level].Nt[tSlice]),makeMatrixOp(data_[level].PatchMWx[slice][np]));
                }


            }

            gsPatchSubassambledLocalOperator<real_t>::Ptr locOperator = gsPatchSubassambledLocalOperator<real_t>::make(localMatricesA,top);
            A[level][slice] = gsParallelOperator<real_t>::make(patchConnection.getConnectionPairs(),locOperator,comm_);

            //     if(output_) gsInfo<<" Finished A!\n";
            if(tSlice>0)
            {
                //      if(output_) gsInfo<<" Start B!\n";
                gsPatchSubassambledLocalOperator<real_t>::Ptr locOperatorB = gsPatchSubassambledLocalOperator<real_t>::make(localMatricesB,top);
                B[level][slice] = gsParallelOperator<real_t>::make(patchConnection.getConnectionPairs(),locOperatorB,comm_);
                //        if(output_) gsInfo<<" Finished B!\n";
            }



        }
    }
    //  if(output_)std::cout<<"finished calc of A,B "<< std::endl;
    //  comm_.barrier();

    ndofs_.resize(data_.back().Mt.size());
    if(myTimeSlices.back().size() == 0)
        gsInfo<<"Error, no Operator on the finest level, something is wrong!\n"<<std::flush;


    //gsInfo<<"Number of dofs on slice: ";
    for(size_t sl = 0; sl< ndofs_.size();++sl)
    {
        //ndofs_[sl] = timeRestr1_.cols()* (spaceLevel == 0? m_prolongation.front()->cols(): m_prolongation[spaceLevel-1]->rows()); //A.back().front()->rows();
        ndofs_[sl] = m_Kx[spaceLevel]->rows()*timeRestr1_.cols();
       // gsInfo<<ndofs_[sl]<<"("<<sl<<")"<<", ";
    }
  //  gsInfo<<"\n";




    real_t maxEigenvalue = -1; //for IETI
    switch(settings_.sliceSolver)
    {
    case Base::FULLSPACETIME:
    {
        precA.resize(TimeLevelsInp_.size());
        LuOfA.resize(TimeLevelsInp_.size());
        for(unsigned int level=0; level<TimeLevelsInp_.size(); level++)
        {
            precA[level].resize(myTimeSlices[TimeLevelsInp_[level]].size());
            LuOfA[level].resize(myTimeSlices[TimeLevelsInp_[level]].size());
            for(size_t slice = 0; slice<myTimeSlices[TimeLevelsInp_[level]].size();++slice)
            {
                int tSlice = myTimeSlices[TimeLevelsInp_[level]][slice];
                //   if(output_) gsInfo<<"start slice "<<tSlice<<"on level "<<level<<"\n";
                switch(settings_.spaceSolver){
                case Base::DIRECT:
                {
                    GISMO_ERROR("Dircet solver not allowed in the space parallel setting");
                    //Note, the solution is distributed if the rhs is distributed!
                    gsMatrix<T> mat;
                    A[level][slice]->toMatrix(mat);

                    precA[level][slice] = makePartialPivLUSolver(mat);
                    LuOfA[level][slice]  = precA[level][slice];
                    break;
                }
                case Base::IETI:
                {

                    /*
                    if(slice==0)
                    {
                        gsMatrix<T> mat;
                        A[level][slice]->toMatrix(mat);
                        gsInfo<<"rank: "<<comm_.rank()<<", Matrix A: \n"<<mat<<"\n\n";

                    }
                    */

                    typename gsIETIAssemblerMPI<T>::Ptr ieti = memory::make_shared(new gsIETIAssemblerMPI<T>(*settings_.m_STAssemblers[spaceLevel],myPatches_,comm_));

                    //    if(output_)std::cout<<"type conversion done for level "<<level<<" and slice "<<slice<<" (mapped to "<<myTimeSlices[TimeLevelsInp_[level]][slice]<<")"<< std::endl;
                    //     comm_.barrier();
		    //settings_.IETIOptions.setString("Scaling","stiff");
                    ieti->setOptions(settings_.IETIOptions);
                    ieti->init();
                    std::vector<gsSparseMatrix<T> > mats(data_[level].PatchMx[slice].size());
                    gsSparseMatrix<T> Mat, temp;
                    for(size_t npi = 0; npi<myPatches_.size();++npi)
                    {
                        size_t np = myPatches_[npi];
                        //        gsInfo<< "first product\n";
                        temp=data_[level].Kt[tSlice]->kron(*data_[level].PatchMx[slice][np]);
                        //         gsInfo<< "second product\n";
                        //Mat=data_[level].Mt[tSlice]->kron(*data_[level].PatchKx[slice][np]);
			// /math::ipow(2,NTimeLevels_-TimeLevelsInp_[level]-1))
			Mat=gsSparseMatrix<T>(*data_[level].Mt[tSlice]).kron(*data_[level].PatchKx[slice][np]);
                        Mat+=temp;
                        //if(slice ==0)
                        // gsInfo<<"rank: "<<comm_.rank()<<", Mat on Patch: "<<np<<"\n "<<Mat.toDense()<<"\n\n";
                        mats[np] = give(Mat);
                    }
                    //     if(output_)std::cout<<"local mats assembled for level "<<level<<" and slice "<<slice<<" (mapped to "<<tSlice<<")"<< std::endl;
                    //     comm_.barrier();
                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                    ieti->assemble();
                    //      if(output_)std::cout<<"IETI setup done for level "<<level<<" and slice "<<slice<<" (mapped to "<<tSlice<<")"<< std::endl<<std::flush;
                    //      comm_.barrier();

                    precA[level][slice]  = IETIAdapterMPI<T>::make(dynamic_cast<gsParallelOperator<T>*>(A[level][slice].get())->getConnectionHandlerTrialSpace(),ieti,m_parHandler,settings_.precondIterations,settings_.precondTol,false,true,maxEigenvalue,outputIETI_iter);
                    //precA[level][slice] = IETIAdapterMPI<T>::make(dynamic_cast<gsParallelOperator<T>*>(A[level][slice].get())->getConnectionHandlerTrialSpace(),ieti,m_parHandler,100,tolExactInnerSolv,false,true,-1,outputIETI_iter);
                    LuOfA[level][slice] = IETIAdapterMPI<T>::make(dynamic_cast<gsParallelOperator<T>*>(A[level][slice].get())->getConnectionHandlerTrialSpace(),ieti,m_parHandler,100,tolExactInnerSolv,false,true,-1,outputIETI_iter);

                   maxEigenvalue = dynamic_cast<IETIAdapterMPI<T>*>(precA[level][slice].get())->getMaxEigenvalue();

                    //     if(output_)std::cout<<"Operator setup done for level "<<level<<" and slice "<<slice<<" (mapped to "<<tSlice<<")"<< std::endl;
                    break;
                }
                };

                //    if(output_)std::cout<<"finished slice "<<tSlice<<"on level "<<level<<"\n";

            }
            //     if(output_)std::cout<<"finished calc of Ainv on level "<<level<< std::endl;
        }
        break;
    case Base::SPATIALWISE:
        {
            precA.resize(TimeLevelsInp_.size());
            LuOfA.resize(TimeLevelsInp_.size());

            for(unsigned int level=0; level<TimeLevelsInp_.size(); level++)
            {
                precA[level].resize(myTimeSlices[TimeLevelsInp_[level]].size());
                LuOfA[level].resize(myTimeSlices[TimeLevelsInp_[level]].size());

                for(size_t slice = 0; slice<myTimeSlices[TimeLevelsInp_[level]].size();++slice)
                {
                    int nx = m_Mx[spaceLevel]->cols();
                    GISMO_ASSERT(m_Mx[spaceLevel]->cols() == m_Mx[spaceLevel]->rows() &&
                                 m_Mx[spaceLevel]->cols() == m_Kx[spaceLevel]->cols() &&
                                 m_Kx[spaceLevel]->cols() == m_Kx[spaceLevel]->rows(), "the spatial matrices have not the same dimension.");

                    int nt = data_[level].Mt[myTimeSlices[TimeLevelsInp_[level]][slice]]->cols();

                    gsInfo<<"on level "<<level<<" and spaceLevel "<<spaceLevel<<" on slice: "<<slice<<" sizes: nx x nt: "<<nx <<" x "<<nt<<"\n";

                    gsMatrix<T> Kt = *data_[level].Kt[myTimeSlices[TimeLevelsInp_[level]][slice]];
                    gsMatrix<T> Mt = *data_[level].Mt[myTimeSlices[TimeLevelsInp_[level]][slice]];

                    switch(settings_.timeDecomp)
                    {
                    case Base::DIAGONALIZATION:
                    {
                        Eigen::GeneralizedEigenSolver<Eigen::Matrix<T,Dynamic,Dynamic> > solv;
                        solv.compute(Kt,Mt,true);
                        //     if(output_)std::cout<<"IETI gave matrices "<<level<<" and slice "<<slice<<" (mapped to "<<tSlice<<")"<< std::endl<<std::flush;
                        //     comm_.barrier();;
                        typename Eigen::GeneralizedEigenSolver<Eigen::Matrix<T,Dynamic,Dynamic> >::EigenvalueType D = solv.eigenvalues();
                        typename gsMatrix<std::complex<T> >::Ptr X = gsMatrix<std::complex<T> > (solv.eigenvectors()).moveToPtr();
                        typename gsMatrix<std::complex<T> >::Ptr VinvTr = gsMatrix<std::complex<T> > ((Mt*(*X))).moveToPtr();
                        //  gsInfo<<"Kt: \n"<<Kt<<"\n\n"<<"Mt: \n"<<Mt<<"\n\n";
                        //   gsInfo<<"X: \n"<<solv.eigenvectors()<<"\n\n"<<"D: \n"<<D<<"\n\n"<<"V^-T: \n"<<*VinvTr<<"\n\n";

                        typename gsBlockOp<std::complex<T> >::Ptr blockPrec = gsBlockOp<std::complex<T> >::make(nt,nt);
                        typename gsBlockOp<std::complex<T> >::Ptr blockEx = gsBlockOp<std::complex<T> >::make(nt,nt);
                        for(int i=0;i<nt;++i)
                        {
                            GISMO_ASSERT(D[i].real()>=0, "The real part of the Eigenvalues is negative, choose a smaller theta");
                            //(K + (alpha+|beta|)M)^{-1}
                            // if(D[i].real()<0)
                            //  gsInfo<<"Warning, EV becomes negative!!!\n\n";
                            T alpha = D[i].real();
                            if(math::abs(D[i].imag()) <complexBound )
                            {
                                switch(settings_.spaceSolver){
                                case Base::DIRECT:
                                {
                                    GISMO_ERROR("Direct Solvers not supported in space-parallel mode");
                                    //typename DistributedComplexify<T>::Ptr solv = DistributedComplexify<T>::make(gsSolverOp<sparseLLTfact>::make((D[i].real()*(*data_[level].Mx[slice]) + *data_[level].Kx[slice])));
                                    //blockPrec->addOperator(i,i,solv);
                                    //blockEx->addOperator(i,i,solv);
                                    break;
                                }
                                case Base::IETI:
                                {
                                    //GISMO_ERROR("IETI not supported for spatialwise decomposition in space-parallel mode");
                                    typename gsIETIAssemblerMPI<T>::Ptr ieti = memory::make_shared(new gsIETIAssemblerMPI<T>(*settings_.m_spatialAssemblers[spaceLevel],myPatches_,comm_));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.resize(data_[level].PatchMx[slice].size());
                                    for(size_t npi = 0; npi<myPatches_.size();++npi)
                                    {
                                        size_t np = myPatches_[npi];
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    alpha * (*data_[level].PatchMx[slice][np])
                                                    + *data_[level].PatchKx[slice][np]);
                                        mats[np]=give(mat);
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();
                                    typename IETIAdapterMPI<T>::Ptr ptr = IETIAdapterMPI<T>::make(m_Kx[spaceLevel].get()->getConnectionHandlerTrialSpace(),ieti,m_parSpaceHandler[spaceLevel],settings_.precondIterations,settings_.precondTol,false,false,maxEigenvalue,outputIETI_iter);
                                    //gsInfo<<"IETI-size: "<<ptr->cols()<<"\n";

                                    blockPrec->addOperator(i,i,gsComplexify<T>::make(ptr));
                                    blockEx->addOperator(i,i,gsComplexify<T>::make(IETIAdapterMPI<T>::make(m_Kx[spaceLevel].get()->getConnectionHandlerTrialSpace(),ieti,m_parSpaceHandler[spaceLevel],100,tolExactInnerSolv,false,false,-1,outputIETI_iter)));
                                    maxEigenvalue = ptr->getMaxEigenvalue();
                                }
                                    break;
                                case Base::MG:
                                {
                                    size_t nLevels = spaceLevel+1;
                                    if(nLevels>1)
                                    {
                                        std::vector<typename gsDistributedOperator<T>::Ptr> linOpsD (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOps (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsR (nLevels-1);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsP (nLevels-1);
                                        for(size_t l=0; l<nLevels;++l)
                                        {
                                            linOpsD[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                            linOps[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                        }
                                        for(size_t l =0; l<nLevels-1;++l)
                                        {
                                            linOpsR[l] = m_restriction[l];
                                            linOpsP[l] = m_prolongation[l];
                                        }

                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);

                                        typename gsMultiGridOp<T>::Ptr mg = gsMultiGridOp<T>::make(linOps,linOpsP,linOpsR,cs);
                                        mg->setOptions(settings_.MGOptions);

                                        for (int l = mg->numLevels() == 1 ? 0 : 1; l < mg->numLevels(); ++l)
                                        {
                                            switch( settings_.smoother ) {
                                            case Smoother_::Richardson:                            mg->setSmoother(l, makeRichardsonOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::Jacobi:                                mg->setSmoother(l, makeJacobiOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::GaussSeidel:                           mg->setSmoother(l, makeGaussSeidelOp(mg->matrix(l))); break;
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD:
                                            {
                                                GISMO_NO_IMPLEMENTATION;
                                                //  std::vector<typename gsLinearOperator<T>::Ptr> localSmoothers = makeSubspaceCorrectedMassSmootherOperatorsDirichlet(settings_.m_spatialAssemblers[i]->multiBasis(),settings_.damping);
                                                break;
                                            }
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:
                                            {
                                                index_t degree = settings_.m_spatialAssemblers[l]->multiBasis().maxCwiseDegree();
                                                const index_t dim = settings_.m_spatialAssemblers[l]->multiBasis().dim();
                                                gsMultiBasis<T> proc_loc_mb;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                    proc_loc_mb.addBasis(settings_.m_spatialAssemblers[l]->multiBasis()[myPatches_[k]].clone().release());

                                                const index_t proc_loc_dofs = (*m_parSpaceHandler[l]).localSize();

                                                gsDofMapper dm; // Global mapper
                                                settings_.m_spatialAssemblers[l]->multiBasis().getMapper(
                                                            (dirichlet::strategy)settings_.m_spatialAssemblers[l]->options().askInt("DirichletStrategy",11),
                                                            (iFace    ::strategy)settings_.m_spatialAssemblers[l]->options().askInt("InterfaceStrategy", 1),
                                                            settings_.m_spatialAssemblers[l]->pde().bc(),
                                                            dm,
                                                            0
                                                            );

                                                // Global => Processor local
                                                gsMatrix<index_t> glob2loc;
                                                {
                                                    gsMatrix<index_t> inp(proc_loc_dofs,1);
                                                    inp.col(0) = gsVector<index_t>::LinSpaced(proc_loc_dofs,0,proc_loc_dofs-1);
                                                    glob2loc.setZero((*m_parSpaceHandler[l]).globalSize(),1);
                                                    (*m_parSpaceHandler[l]).addLocalVectorToGlobal( inp, glob2loc );
                                                }

                                                std::vector< gsVector<index_t> > proc_loc_maps;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                {
                                                    const index_t sz = proc_loc_mb[k].size();
                                                    index_t kk = myPatches_[k];
                                                    gsVector<index_t> local(sz);
                                                    for (index_t j=0; j<sz; ++j)
                                                    {
                                                        local[j] = dm.is_free(j,kk) ? glob2loc(dm.index(j, kk),0) : -1u;
                                                    }
                                                    proc_loc_maps.push_back(give(local));
                                                }
                                                std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<> > > > pieces =
                                                        constructPieces( proc_loc_mb, proc_loc_maps, proc_loc_dofs );

                                                std::vector< typename gsLinearOperator<T>::Ptr > localSmoothers;
                                                std::vector< gsSparseMatrix<T> > smootherTransfers;
                                                const index_t nrPieces = pieces.size();

                                                real_t h = 1;
                                                for ( size_t jj=0; jj<settings_.m_spatialAssemblers[l]->multiBasis().nBases(); ++jj)
                                                {
                                                    const gsGeometry<T>& patch = settings_.m_spatialAssemblers[l]->patches().patch(jj);
                                                    real_t diam =(patch.coefAtCorner(boxCorner::getFirst(dim)) - patch.coefAtCorner(boxCorner::getLast(dim))).norm();
                                                    h = std::min(h,settings_.m_spatialAssemblers[l]->multiBasis()[jj].getMinCellLength()*diam);
                                                }

                                                for ( index_t dd = 0; dd<nrPieces; ++dd )
                                                {
                                                    gsBoundaryConditions<T> localbc;
                                                    for( index_t ps=0; ps < 2*dim; ++ps )
                                                        localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

                                                    const index_t sz = pieces[dd].size();
                                                    for ( index_t jj=0; jj<sz; ++jj)
                                                    {
                                                        if (dd == dim)
                                                        {
                                                            localSmoothers.push_back( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc,alpha) );
                                                        }
                                                        else if (dd > 1)
                                                        {
                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            localSmoothers.push_back( gsScaledOp<T>::make( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc,(scalingM+alpha*math::pow(2,dim-dd)*math::pow(h/(2*degree+1),dim-dd) )/scalingK), 1/scalingK ) );
                                                        }
                                                        else if (dd == 1)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );

                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            gsSparseMatrix<T> M, K;
                                                            assembleParameterMass(*pieces[dd][jj].first, M);
                                                            assembleParameterStiffness(*pieces[dd][jj].first, K);
                                                            M = M.block( 1, 1, M.rows()-2, M.cols()-2 );
                                                            K = K.block( 1, 1, K.rows()-2, K.cols()-2 );
                                                            gsSparseMatrix<T> mat = scalingK * K + (scalingM+alpha*math::pow(2,dim-dd)*math::pow( h/(2*degree+1), dim-dd )) * M;

                                                            localSmoothers.push_back( makeSparseCholeskySolver(mat) );
                                                        }
                                                        else if (dd == 0)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );
                                                            const index_t sz = pieces[dd][jj].second.rows();
                                                            const real_t scalingFactor = math::pow(2,dim) * dim * math::pow( h/(1+2*degree), dim-1 ) * degree*degree/(h*(2*degree-1))
                                                                    + math::pow(h/(2*degree+1),dim)*math::pow(2,dim)*alpha;
                                                            localSmoothers.push_back( gsScaledOp<T>::make( gsIdentityOp<T>::make(sz), 1/scalingFactor) );
                                                        }

                                                        smootherTransfers.push_back( give( pieces[dd][jj].second ) );
                                                    }
                                                }

                                                mg->setSmoother(
                                                            l,
                                                            gsParallelAdditivePreconditionerOp<T>::make(
                                                                linOpsD[l],
                                                                give(smootherTransfers),
                                                                give(localSmoothers),
                                                                settings_.outerDamping
                                                                )
                                                            );

                                                break;
                                            }


                                            }
                                        }

                                        typename gsParallelPreconditionerOp<T>::Ptr mg_prec = gsParallelPreconditionerOp<T>::make(mg);

                                        blockPrec->addOperator(i,i,gsComplexify<T>::make(mg_prec));
                                        gsOptionList options;
                                        options.addReal("Iterations","number of maximal Iterations", 100);
                                        options.addReal("Tolerance", "tolerance", tolExactInnerSolv);
                                        blockEx->addOperator(i,i,gsComplexify<T>::make(gsDistributedIterativeSolverOp<gsParallelCG<T> >::make(m_Kx.back()->getConnectionHandlerTestSpace(),linOpsD.back(), mg_prec,options)));
                                    }
                                    else
                                    {
                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);
                                        blockPrec->addOperator(i,i,gsComplexify<T>::make(cs));
                                        blockEx->addOperator(i,i,gsComplexify<T>::make(cs));

                                    }
                                }
                                    break;
                                };
                            }
                            else
                            {
                                alpha = D[i].real()+math::abs(D[i].imag());
                                for(int j=i+1;j<nt;++j)
                                {
                                    if(std::abs(D[i] - std::conj(D[j]))<2.e-10)
                                    {
                                        typename gsDistributedOperator<T>::Ptr precond;

                                        switch(settings_.spaceSolver){
                                        case Base::DIRECT:
                                            GISMO_ERROR("Direct Solvers not supported in space-parallel mode");
                                            //precond= gsSolverOp<sparseLLTfact>::make((gsSparseMatrix<T>(*data_[level].Kx[slice]+(D[i].real()+math::abs(D[i].imag()))**data_[level].Mx[slice])));
                                            break;
                                        case Base::IETI:
                                        {
                                            typename gsIETIAssemblerMPI<T>::Ptr ieti = memory::make_shared(new gsIETIAssemblerMPI<T>(*settings_.m_spatialAssemblers[spaceLevel],myPatches_,comm_));
                                            ieti->setOptions(settings_.IETIOptions);
                                            ieti->init();
                                            std::vector<gsSparseMatrix<T> > mats;
                                            mats.resize(data_[level].PatchMx[slice].size());
                                            for(size_t npi = 0; npi<myPatches_.size();++npi)
                                            {
                                                size_t np = myPatches_[npi];
                                                gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                            alpha* (*data_[level].PatchMx[slice][np]) + *data_[level].PatchKx[slice][np]);
                                                mats[np]=give(mat);
                                            }
                                            ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                            ieti->assemble();
                                            precond = IETIAdapterMPI<T>::make(m_Kx[spaceLevel].get()->getConnectionHandlerTrialSpace(),ieti,m_parSpaceHandler[spaceLevel],settings_.precondIterations,settings_.precondTol,false,false,maxEigenvalue,outputIETI_iter);
                                            maxEigenvalue =dynamic_cast<IETIAdapterMPI<T>*>(precond.get())->getMaxEigenvalue();
                                        }
                                            break;
                                        case Base::MG:
                                        {
                                            size_t nLevels = spaceLevel+1;
                                            if(nLevels>1)
                                            {
                                                std::vector<typename gsDistributedOperator<T>::Ptr> linOpsD (nLevels);
                                                std::vector<typename gsLinearOperator<T>::Ptr> linOps (nLevels);
                                                std::vector<typename gsLinearOperator<T>::Ptr> linOpsR (nLevels-1);
                                                std::vector<typename gsLinearOperator<T>::Ptr> linOpsP (nLevels-1);
                                                for(size_t l=0; l<nLevels;++l)
                                                {
                                                    linOpsD[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                                    linOps[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                                }
                                                for(size_t l =0; l<nLevels-1;++l)
                                                {
                                                    linOpsR[l] = m_restriction[l];
                                                    linOpsP[l] = m_prolongation[l];
                                                }

                                                std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                                for(size_t k=0; k<myPatches_.size();++k)
                                                {
                                                    size_t np = myPatches_[k];
                                                    c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                                }

                                                typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);

                                                typename gsMultiGridOp<T>::Ptr mg = gsMultiGridOp<T>::make(linOps,linOpsP,linOpsR,cs );
                                                mg->setOptions(settings_.MGOptions);

                                                for (int l = mg->numLevels() == 1 ? 0 : 1; l < mg->numLevels(); ++l)
                                                {
                                                    switch( settings_.smoother ) {
                                                    case Smoother_::Richardson:                            mg->setSmoother(l, makeRichardsonOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                                    case Smoother_::Jacobi:                                mg->setSmoother(l, makeJacobiOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                                    case Smoother_::GaussSeidel:                           mg->setSmoother(l, makeGaussSeidelOp(mg->matrix(l))); break;
                                                    case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD:
                                                    {
                                                        std::vector<typename gsLinearOperator<T>::Ptr> localSmoothers = makeSubspaceCorrectedMassSmootherOperatorsDirichlet(settings_.m_spatialAssemblers[l]->multiBasis(),settings_.damping);
                                                        GISMO_NO_IMPLEMENTATION;
                                                        /*
                                                     mg->setSmoother(
                                                                 i,
                                                                 gsParallelAdditivePreconditionerOp<T>::make(
                                                                     linOpsD[i],
                                                                     setupPiecewisePreconditioner<T>(
                                                                         *linOpsD[i], //NeedsParallelOperator!!!
                                                                         *m_parSpaceHandler[i],
                                                                         give(localSmoothers),
                                                                         settings_.m_spatialAssemblers[i]->patches(),
                                                                         settings_.m_spatialAssemblers[i]->multiBasis(),
                                                                         settings_.m_spatialAssemblers[i]->pde().bc(),
                                                                         settings_.m_spatialAssemblers[i]->options(),
                                                                         myPatches_,
                                                                         comm_
                                                                         ),
                                                                     settings_.outerDamping
                                                                     )
                                                                 );
                                                                 */
                                                        break;
                                                    }
                                                    case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:
                                                    {
                                                        index_t degree = settings_.m_spatialAssemblers[l]->multiBasis().maxCwiseDegree();
                                                        const index_t dim = settings_.m_spatialAssemblers[l]->multiBasis().dim();
                                                        gsMultiBasis<T> proc_loc_mb;
                                                        for (size_t k=0;k<myPatches_.size();++k)
                                                            proc_loc_mb.addBasis(settings_.m_spatialAssemblers[l]->multiBasis()[myPatches_[k]].clone().release());

                                                        const index_t proc_loc_dofs = (*m_parSpaceHandler[l]).localSize();

                                                        gsDofMapper dm; // Global mapper
                                                        settings_.m_spatialAssemblers[l]->multiBasis().getMapper(
                                                                    (dirichlet::strategy)settings_.m_spatialAssemblers[l]->options().askInt("DirichletStrategy",11),
                                                                    (iFace    ::strategy)settings_.m_spatialAssemblers[l]->options().askInt("InterfaceStrategy", 1),
                                                                    settings_.m_spatialAssemblers[l]->pde().bc(),
                                                                    dm,
                                                                    0
                                                                    );

                                                        // Global => Processor local
                                                        gsMatrix<index_t> glob2loc;
                                                        {
                                                            gsMatrix<index_t> inp(proc_loc_dofs,1);
                                                            inp.col(0) = gsVector<index_t>::LinSpaced(proc_loc_dofs,0,proc_loc_dofs-1);
                                                            glob2loc.setZero((*m_parSpaceHandler[l]).globalSize(),1);
                                                            (*m_parSpaceHandler[l]).addLocalVectorToGlobal( inp, glob2loc );
                                                        }

                                                        std::vector< gsVector<index_t> > proc_loc_maps;
                                                        for (size_t k=0;k<myPatches_.size();++k)
                                                        {
                                                            const index_t sz = proc_loc_mb[k].size();
                                                            index_t kk = myPatches_[k];
                                                            gsVector<index_t> local(sz);
                                                            for (index_t jj=0; jj<sz; ++jj)
                                                            {
                                                                local[jj] = dm.is_free(jj,kk) ? glob2loc(dm.index(jj, kk),0) : -1u;
                                                            }
                                                            proc_loc_maps.push_back(give(local));
                                                        }
                                                        std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<> > > > pieces =
                                                                constructPieces( proc_loc_mb, proc_loc_maps, proc_loc_dofs );

                                                        std::vector< typename gsLinearOperator<T>::Ptr > localSmoothers;
                                                        std::vector< gsSparseMatrix<T> > smootherTransfers;
                                                        const index_t nrPieces = pieces.size();

                                                        real_t h = 1;
                                                        for ( size_t jj=0; jj<settings_.m_spatialAssemblers[l]->multiBasis().nBases(); ++jj)
                                                        {
                                                            const gsGeometry<T>& patch = settings_.m_spatialAssemblers[l]->patches().patch(jj);
                                                            real_t diam =(patch.coefAtCorner(boxCorner::getFirst(dim)) - patch.coefAtCorner(boxCorner::getLast(dim))).norm();
                                                            h = std::min(h,settings_.m_spatialAssemblers[l]->multiBasis()[jj].getMinCellLength()*diam);
                                                        }
                                                        for ( index_t dd = 0; dd<nrPieces; ++dd )
                                                        {
                                                            gsBoundaryConditions<T> localbc;
                                                            for( index_t ps=0; ps < 2*dim; ++ps )
                                                                localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

                                                            const index_t sz = pieces[dd].size();
                                                            for ( index_t jj=0; jj<sz; ++jj)
                                                            {
                                                                if (dd == dim)
                                                                {
                                                                    localSmoothers.push_back( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc,alpha));
                                                                }
                                                                else if (dd > 1)
                                                                {
                                                                    const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                                    const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                                    localSmoothers.push_back( gsScaledOp<T>::make( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc, (scalingM+alpha*math::pow(2,dim-dd)*math::pow(h/(2*degree+1),dim-dd) )/scalingK), 1/scalingK ) );
                                                                }
                                                                else if (dd == 1)
                                                                {
                                                                    //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                                    //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );

                                                                    const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                                    const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                                    gsSparseMatrix<T> M, K;
                                                                    assembleParameterMass(*pieces[dd][jj].first, M);
                                                                    assembleParameterStiffness(*pieces[dd][jj].first, K);
                                                                    M = M.block( 1, 1, M.rows()-2, M.cols()-2 );
                                                                    K = K.block( 1, 1, K.rows()-2, K.cols()-2 );
                                                                    gsSparseMatrix<T> mat = scalingK * K + (scalingM+ alpha*math::pow(2,dim-dd)*math::pow( h/(2*degree+1), dim-dd )) * M;

                                                                    localSmoothers.push_back( makeSparseCholeskySolver(mat) );
                                                                }
                                                                else if (dd == 0)
                                                                {
                                                                    //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                                    //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );
                                                                    const index_t sz = pieces[dd][jj].second.rows();
                                                                    const real_t scalingFactor = math::pow(2,dim) * dim * math::pow( h/(1+2*degree), dim-1 ) * degree*degree/(h*(2*degree-1))+ math::pow(h/(2*degree+1),dim)*math::pow(2,dim)*alpha;
                                                                    localSmoothers.push_back( gsScaledOp<T>::make( gsIdentityOp<T>::make(sz), 1/scalingFactor) );
                                                                }

                                                                smootherTransfers.push_back( give( pieces[dd][jj].second ) );
                                                            }
                                                        }

                                                        mg->setSmoother(
                                                                    l,
                                                                    gsParallelAdditivePreconditionerOp<T>::make(
                                                                        linOpsD[l],
                                                                        give(smootherTransfers),
                                                                        give(localSmoothers),
                                                                        settings_.outerDamping
                                                                        )
                                                                    );
                                                        break;
                                                    }


                                                    }
                                                }

                                                precond = gsParallelPreconditionerOp<T>::make(mg);
                                               // gsInfo<<"precond size: "<<precond->rows()<<"  -  "<<"mg size: "<<mg->rows()<<"\n";
                                            }
                                            else
                                            {
                                              //  gsInfo<<"Calling coarse grid solver only:\n\n";
                                                std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                                for(size_t k=0; k<myPatches_.size();++k)
                                                {
                                                    size_t np = myPatches_[k];
                                                    c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                                }

                                                typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);
                                                precond = gsParallelOperator<T>::make(m_Kx[spaceLevel]->getConnectionHandlerTrialSpace(),m_Kx[spaceLevel]->getConnectionHandlerTestSpace(),cs);
                                              //  gsInfo<<"precond size: "<<precond->rows()<<"  -  "<<"mg size: "<<cs->rows()<<"\n";
                                            }

                                            break;
                                        }
                                        };

                                        /*
                                        typename gsMatrix<T>::EigenSolver a;
                                        gsMatrix<T> locMat,locMat_, Prec,col,col_;
                                        gsMatrix<T> K,M;
                                        gsMatrix<T> MK;
                                        if(slice == 0)
                                        {
                                            gsMatrix<T> id = gsMatrix<T>::Identity(m_parSpaceHandler[spaceLevel]->globalSize(),m_parSpaceHandler[spaceLevel]->globalSize());
                                            m_parSpaceHandler[spaceLevel]->extractLocalVector(id,locMat_);
                                            precond->distribute(locMat_,locMat);
                                            for(int c=0; c<locMat.cols();++c)
                                            {
                                                precond->apply(locMat.col(c),col);
                                                precond->distribute(col,col_);
                                                locMat.col(c) = col_;
                                            }
                                            m_parSpaceHandler[spaceLevel]->buildGlobalVector(locMat,Prec);
                                            a.compute(Prec);
                                            T symmNorm = (Prec-Prec.transpose()).norm();
                                            if(comm_.rank()==0) gsInfo<<"\t\t MG: eigenvalues for lambda: "<<D[i]<<" -> "<<a.eigenvalues().real().minCoeff()<<" - "<<a.eigenvalues().real().maxCoeff()<<" || symm: "<<symmNorm<<" \n";


                                            K= m_Kx[spaceLevel]->toDenseMatrix(*m_parSpaceHandler[spaceLevel]);
                                            M= m_Mx[spaceLevel]->toDenseMatrix(*m_parSpaceHandler[spaceLevel]);

                                            MK = K + alpha*M;
                                            symmNorm = (MK-MK.transpose()).norm();
                                            a.compute(MK);
                                            if(comm_.rank()==0) gsInfo<<"\t\t KM eigenvalues for lambda: "<<D[i]<<" -> "<<a.eigenvalues().real().minCoeff()<<" - "<<a.eigenvalues().real().maxCoeff()<<" || symm: "<<symmNorm<<" \n";
                                            a.compute(Prec*MK);
                                            if(comm_.rank()==0) gsInfo<<"\t\t MG-KM: eigenvalues for lambda: "<<D[i]<<" -> "<<a.eigenvalues().real().minCoeff()<<" - "<<a.eigenvalues().real().maxCoeff()<<"\n"<<std::flush;
                                        }
*/
                                        /*
                                        for(int z=0; z<30;z++)
                                        {
                                            typename gsDistributedOperator<T>::Ptr op = gsDistributedSumOp<T>::make(m_Kx[spaceLevel],gsDistributedScaledOp<T>::make(m_Mx[spaceLevel],alpha));
                                            gsParallelCG<T> cg(op,precond,comm_);
                                            cg.setMaxIterations(200);
                                            cg.setCalcEigenvalues(true);
                                            cg.setTolerance(1.e-8);
                                            gsMatrix<T> rhs = gsMatrix<T>::Random(op->rows(),1);
                                            gsMatrix<T> x = gsMatrix<T>::Zero(op->rows(),1);
                                            gsMatrix<T> hist;
                                            cg.solveDetailed(rhs,x,hist);
                                            gsMatrix<T> eigs;
                                            cg.getEigenvalues(eigs);
                                           // gsInfo<<"Hist: \n "<<hist<<"\n";
                                            gsInfo<<"cg for alpha+|beta| converged after: "<<cg.iterations()<<" iterations with residual: "<<cg.error()<<" and ev: "<<eigs.minCoeff()<<" - "<<eigs.maxCoeff()<<"\n"<<std::flush;;
                                        }
                                        gsInfo<<"\n"<<std::flush;
*/


                                        blockPrec->addOperator(i,i,ComplexSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],D[i],precond,comm_,settings_.approxIterates,settings_.tol,true));
                                        blockEx->addOperator(i,i,ComplexSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],D[i],precond,comm_,100,tolExactInnerSolv,true)); //TODO: false does not work here

                                        blockPrec->addOperator(j,j,ComplexSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],D[j],precond,comm_,settings_.approxIterates,settings_.tol,true));
                                        blockEx->addOperator(j,j,ComplexSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],D[j],precond,comm_,100,tolExactInnerSolv,true));//TODO: false does not work here
                                        break;
                                    }
                                }
                            }
                        }

                        typename gsLinearOperator<std::complex<T> >::Ptr V = makePartialPivLUSolver(VinvTr);
                      //  gsInfo<<"blockPrec size: "<<blockPrec->rows()<<" x "<<blockPrec->cols()<<" --- Kron: "<<V->rows()*nx<<" x "<<V->cols()*nx<<"\n";


                        precA[level][slice] = gsDeComplexify<T>::make( gsProductOp<std::complex<T> >::make(gsKroneckerOp<std::complex<T> >::make(V,gsIdentityOp<std::complex<T> >::make(nx)),blockPrec,gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(X),gsIdentityOp<std::complex<T> >::make(nx))));
                        LuOfA[level][slice] = gsDeComplexify<T>::make( gsProductOp<std::complex<T> >::make(gsKroneckerOp<std::complex<T> >::make(V,gsIdentityOp<std::complex<T> >::make(nx)),blockEx,gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(X),gsIdentityOp<std::complex<T> >::make(nx))));

                    }
                        break;
                    case Base::COMPLEXSCHUR:
                    {
                        gsMatrix<T> Minv = Mt.inverse();
                        Eigen::ComplexSchur<Eigen::Matrix<T,Dynamic,Dynamic> > solv(Minv*Kt, true);
                        typename gsMatrix<std::complex<T> >::Ptr Pinv =  (gsMatrix<std::complex<T> >(solv.matrixU().adjoint()*Minv)).moveToPtr();
                        typename gsMatrix<std::complex<T> >::Ptr Q = (gsMatrix<std::complex<T> >(solv.matrixU())).moveToPtr();
                        typename gsMatrix<std::complex<T> >::Ptr UT = (gsMatrix<std::complex<T> >(solv.matrixT())).moveToPtr();

                        std::vector<typename gsLinearOperator<std::complex<T> >::Ptr> APrec(nt);
                        std::vector<typename gsLinearOperator<std::complex<T> >::Ptr> AEx(nt);


                        for(int i=0;i<nt;++i)
                        {
                            GISMO_ASSERT((*UT)(i,i).real()>=0, "The real part of the Eigenvalues is negative, choose a smaller theta");
                            //(K + (alpha+|beta|)M)^{-1}
                            if(math::abs((*UT)(i,i).imag()) <complexBound)
                            {
                                switch(settings_.spaceSolver){
                                case Base::DIRECT:
                                    GISMO_ERROR("Direct Solvers not supported in space-parallel mode");
                                    //APrec[i]=Complexify<T>::make(gsSolverOp<sparseLLTfact>::make(((*UT)(i,i).real()*(*data_[level].Mx[slice]) + *data_[level].Kx[slice])));
                                    //AEx[i]=APrec[i];
                                    break;
                                case Base::IETI:
                                {
                                    typename gsIETIAssemblerMPI<T>::Ptr ieti = memory::make_shared(new gsIETIAssemblerMPI<T>(*settings_.m_spatialAssemblers[spaceLevel],myPatches_,comm_));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.reserve(data_[level].PatchMx[slice].size());
                                    for(size_t npi = 0; npi<myPatches_.size();++npi)
                                    {
                                        size_t np = myPatches_[npi];
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    (*UT)(i, i).real() * (*data_[level].PatchMx[slice][np])
                                                    + *data_[level].PatchKx[slice][np]);
                                        mats.push_back(give(mat));
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();

                                    typename IETIAdapterMPI<T>::Ptr ptr = IETIAdapterMPI<T>::make(m_Kx[spaceLevel].get()->getConnectionHandlerTrialSpace(),ieti,m_parSpaceHandler[spaceLevel],settings_.precondIterations,settings_.precondTol,false,false,-1,outputIETI_iter);
                                    //gsInfo<<"IETI-size: "<<ptr->cols()<<"\n";

                                    APrec[i]=gsComplexify<T>::make(ptr);
                                    AEx[i]  = gsComplexify<T>::make(IETIAdapterMPI<T>::make(m_Kx[spaceLevel].get()->getConnectionHandlerTrialSpace(),ieti,m_parSpaceHandler[spaceLevel],100,tolExactInnerSolv,false,false,-1,outputIETI_iter));
                                    //maxEigenvalue = ptr->getMaxEigenvalue();
                                }
                                    break;
                                case Base::MG:
                                {
                                    size_t nLevels = spaceLevel+1;
                                    if(nLevels>1)
                                    {
                                        std::vector<typename gsDistributedOperator<T>::Ptr> linOpsD (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOps (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsR (nLevels-1);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsP (nLevels-1);
                                        for(size_t l=0; l<nLevels;++l)
                                        {
                                            linOpsD[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],(*UT)(i, i).real()));
                                            linOps[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],(*UT)(i, i).real()));
                                        }
                                        for(size_t l =0; l<nLevels-1;++l)
                                        {
                                            linOpsR[l] = m_restriction[l];
                                            linOpsP[l] = m_prolongation[l];
                                        }

                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*(*UT)(i, i).real()).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);

                                        typename gsMultiGridOp<T>::Ptr mg = gsMultiGridOp<T>::make(linOps,linOpsP,linOpsR,cs);
                                        mg->setOptions(settings_.MGOptions);

                                        for (int l = mg->numLevels() == 1 ? 0 : 1; l < mg->numLevels(); ++l)
                                        {
                                            switch( settings_.smoother ) {
                                            case Smoother_::Richardson:                            mg->setSmoother(l, makeRichardsonOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::Jacobi:                                mg->setSmoother(l, makeJacobiOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::GaussSeidel:                           mg->setSmoother(l, makeGaussSeidelOp(mg->matrix(l))); break;
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD:
                                            {
                                                GISMO_NO_IMPLEMENTATION;
                                                //  std::vector<typename gsLinearOperator<T>::Ptr> localSmoothers = makeSubspaceCorrectedMassSmootherOperatorsDirichlet(settings_.m_spatialAssemblers[i]->multiBasis(),settings_.damping);
                                                /*
                                             mg->setSmoother(
                                                         i,
                                                         gsParallelAdditivePreconditionerOp<T>::make(
                                                             linOpsD[i],
                                                             setupPiecewisePreconditioner<T>(
                                                                 *linOpsD[i], //NeedsParallelOperator!!!
                                                                 *m_parSpaceHandler[i],
                                                                 give(localSmoothers),
                                                                 settings_.m_spatialAssemblers[i]->patches(),
                                                                 settings_.m_spatialAssemblers[i]->multiBasis(),
                                                                 settings_.m_spatialAssemblers[i]->pde().bc(),
                                                                 settings_.m_spatialAssemblers[i]->options(),
                                                                 myPatches_,
                                                                 comm_
                                                                 ),
                                                             settings_.outerDamping
                                                             )
                                                         );
                                                         */
                                                break;
                                            }
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:
                                            {
                                                index_t degree = settings_.m_spatialAssemblers[l]->multiBasis().maxCwiseDegree();
                                                const index_t dim = settings_.m_spatialAssemblers[l]->multiBasis().dim();
                                                gsMultiBasis<T> proc_loc_mb;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                    proc_loc_mb.addBasis(settings_.m_spatialAssemblers[l]->multiBasis()[myPatches_[k]].clone().release());

                                                const index_t proc_loc_dofs = (*m_parSpaceHandler[l]).localSize();

                                                gsDofMapper dm; // Global mapper
                                                settings_.m_spatialAssemblers[l]->multiBasis().getMapper(
                                                            (dirichlet::strategy)settings_.m_spatialAssemblers[l]->options().askInt("DirichletStrategy",11),
                                                            (iFace    ::strategy)settings_.m_spatialAssemblers[l]->options().askInt("InterfaceStrategy", 1),
                                                            settings_.m_spatialAssemblers[l]->pde().bc(),
                                                            dm,
                                                            0
                                                            );

                                                // Global => Processor local
                                                gsMatrix<index_t> glob2loc;
                                                {
                                                    gsMatrix<index_t> inp(proc_loc_dofs,1);
                                                    inp.col(0) = gsVector<index_t>::LinSpaced(proc_loc_dofs,0,proc_loc_dofs-1);
                                                    glob2loc.setZero((*m_parSpaceHandler[l]).globalSize(),1);
                                                    (*m_parSpaceHandler[l]).addLocalVectorToGlobal( inp, glob2loc );
                                                }

                                                std::vector< gsVector<index_t> > proc_loc_maps;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                {
                                                    const index_t sz = proc_loc_mb[k].size();
                                                    index_t kk = myPatches_[k];
                                                    gsVector<index_t> local(sz);
                                                    for (index_t j=0; j<sz; ++j)
                                                    {
                                                        local[j] = dm.is_free(j,kk) ? glob2loc(dm.index(j, kk),0) : -1u;
                                                    }
                                                    proc_loc_maps.push_back(give(local));
                                                }
                                                std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<> > > > pieces =
                                                        constructPieces( proc_loc_mb, proc_loc_maps, proc_loc_dofs );

                                                std::vector< typename gsLinearOperator<T>::Ptr > localSmoothers;
                                                std::vector< gsSparseMatrix<T> > smootherTransfers;
                                                const index_t nrPieces = pieces.size();

                                                real_t h = 1;
                                                for ( size_t jj=0; jj<settings_.m_spatialAssemblers[l]->multiBasis().nBases(); ++jj)
                                                {
                                                    const gsGeometry<T>& patch = settings_.m_spatialAssemblers[l]->patches().patch(jj);
                                                    real_t diam =(patch.coefAtCorner(boxCorner::getFirst(dim)) - patch.coefAtCorner(boxCorner::getLast(dim))).norm();
                                                    h = std::min(h,settings_.m_spatialAssemblers[l]->multiBasis()[jj].getMinCellLength()*diam);
                                                }
                                                for ( index_t dd = 0; dd<nrPieces; ++dd )
                                                {
                                                    gsBoundaryConditions<T> localbc;
                                                    for( index_t ps=0; ps < 2*dim; ++ps )
                                                        localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

                                                    const index_t sz = pieces[dd].size();
                                                    for ( index_t jj=0; jj<sz; ++jj)
                                                    {
                                                        if (dd == dim)
                                                        {
                                                            localSmoothers.push_back( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc,(*UT)(i,i).real()) );
                                                        }
                                                        else if (dd > 1)
                                                        {
                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            localSmoothers.push_back( gsScaledOp<T>::make( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc, (scalingM+(*UT)(i,i).real()*math::pow(2,dim-dd)*math::pow(h/(2*degree+1),dim-dd) )/scalingK), 1/scalingK ) );
                                                        }
                                                        else if (dd == 1)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );

                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            gsSparseMatrix<T> M, K;
                                                            assembleParameterMass(*pieces[dd][jj].first, M);
                                                            assembleParameterStiffness(*pieces[dd][jj].first, K);
                                                            M = M.block( 1, 1, M.rows()-2, M.cols()-2 );
                                                            K = K.block( 1, 1, K.rows()-2, K.cols()-2 );
                                                            gsSparseMatrix<T> mat = scalingK * K + (scalingM+(*UT)(i,i).real()*math::pow(2,dim-dd)*math::pow( h/(2*degree+1), dim-dd )) * M;

                                                            localSmoothers.push_back( makeSparseCholeskySolver(mat) );
                                                        }
                                                        else if (dd == 0)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );
                                                            const index_t sz = pieces[dd][jj].second.rows();
                                                            const real_t scalingFactor = math::pow(2,dim) * dim * math::pow( h/(1+2*degree), dim-1 ) * degree*degree/(h*(2*degree-1))
                                                                    + math::pow(h/(2*degree+1),dim)*math::pow(2,dim)*(*UT)(i,i).real();
                                                            localSmoothers.push_back( gsScaledOp<T>::make( gsIdentityOp<T>::make(sz), 1/scalingFactor) );
                                                        }

                                                        smootherTransfers.push_back( give( pieces[dd][jj].second ) );
                                                    }
                                                }

                                                mg->setSmoother(
                                                            l,
                                                            gsParallelAdditivePreconditionerOp<T>::make(
                                                                linOpsD[l],
                                                                give(smootherTransfers),
                                                                give(localSmoothers),
                                                                settings_.outerDamping
                                                                )
                                                            );
                                                break;
                                            }
                                            }
                                        }

                                        typename gsParallelPreconditionerOp<T>::Ptr mg_prec = gsParallelPreconditionerOp<T>::make(mg);

                                        APrec[i]  = gsComplexify<T>::make(mg_prec);
                                        gsOptionList options;
                                        options.addReal("Iterations","number of maximal Iterations", 100);
                                        options.addReal("Tolerance", "tolerance", tolExactInnerSolv);
                                        AEx[i] = gsComplexify<T>::make(gsDistributedIterativeSolverOp<gsParallelCG<T> >::make(m_Kx.back()->getConnectionHandlerTestSpace(),linOpsD.back(), mg_prec,options));

                                        /*
                                        typename gsMatrix<T>::EigenSolver a;
                                        gsMatrix<T> locMat,locMat_, Prec,col,col_;
                                        gsMatrix<T> K,M;
                                        gsMatrix<T> MK;
                                        if(slice == 0)
                                        {
                                            gsMatrix<T> id = gsMatrix<T>::Identity(m_parSpaceHandler[spaceLevel]->globalSize(),m_parSpaceHandler[spaceLevel]->globalSize());

                                            m_parSpaceHandler[spaceLevel]->extractLocalVector(id,locMat_);
                                            mg_prec->distribute(locMat_,locMat);
                                            for(int c=0; c<locMat.cols();++c)
                                            {
                                                mg_prec->apply(locMat.col(c),col_);
                                                mg_prec->distribute(col_,col);
                                                locMat.col(c) = col;
                                            }

                                            m_parSpaceHandler[spaceLevel]->buildGlobalVector(locMat,Prec);
                                            T symmNorm = (Prec-Prec.transpose()).norm();
                                            a.compute(Prec);
                                            if(comm_.rank()==0) gsInfo<<"\t\t MG: eigenvalues for lambda: "<<(*UT)(i,i)<<" -> "<<a.eigenvalues().real().minCoeff()<<" - "<<a.eigenvalues().real().maxCoeff()<<" || symm: "<<symmNorm<<" \n";


                                            K= m_Kx[spaceLevel]->toDenseMatrix(*m_parSpaceHandler[spaceLevel]);
                                            M= m_Mx[spaceLevel]->toDenseMatrix(*m_parSpaceHandler[spaceLevel]);

                                            MK = K + ((*UT)(i,i)).real()*M;
                                            symmNorm = (MK-MK.transpose()).norm();
                                            a.compute(MK);
                                            if(comm_.rank()==0) gsInfo<<"\t\t KM eigenvalues for lambda: "<<(*UT)(i,i)<<" -> "<<a.eigenvalues().real().minCoeff()<<" - "<<a.eigenvalues().real().maxCoeff()<<" || symm: "<<symmNorm<<" \n";
                                            a.compute(Prec*MK);
                                            symmNorm = (Prec*MK-(Prec*MK).transpose()).norm();
                                            if(comm_.rank()==0) gsInfo<<"\t\t MG-KM: eigenvalues for lambda: "<<(*UT)(i,i)<<" -> "<<a.eigenvalues().real().minCoeff()<<" - "<<a.eigenvalues().real().maxCoeff()<<" || symm: "<<symmNorm<<" \n";;
                                        }
                                        */
                                    }
                                    else
                                    {
                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*(*UT)(i, i).real()).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);

                                        APrec[i]  =  gsComplexify<T>::make(cs);
                                        AEx[i] = gsComplexify<T>::make(cs);

                                    }
                                }
                                    break;
                                default:
                                    GISMO_NO_IMPLEMENTATION;
                                };
                            }
                            else if(APrec[i]==NULL)
                            {
                                T alpha = (*UT)(i, i).real() + math::abs((*UT)(i, i).imag());
                                typename gsDistributedOperator<T>::Ptr precond;

                                switch(settings_.spaceSolver){
                                case Base::DIRECT:
                                    GISMO_ERROR("Direct Solvers not supported in space-parallel mode");
                                    //precond= gsSolverOp<sparseLLTfact>::make((gsSparseMatrix<T>(*data_[level].Kx[slice]+((*UT)(i,i).real()+math::abs((*UT)(i,i).imag()))**data_[level].Mx[slice])));
                                    break;
                                case Base::IETI:
                                {
                                    typename gsIETIAssemblerMPI<T>::Ptr ieti = memory::make_shared(new gsIETIAssemblerMPI<T>(*settings_.m_spatialAssemblers[spaceLevel],myPatches_,comm_));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.reserve(data_[level].PatchMx[slice].size());
                                    for(size_t npi = 0; npi<myPatches_.size();++npi)
                                    {
                                        size_t np = myPatches_[npi];
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    alpha * (*data_[level].PatchMx[slice][np]) + *data_[level].PatchKx[slice][np]);
                                        mats.push_back(give(mat));
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();

                                    precond = IETIAdapterMPI<T>::make(m_Kx[spaceLevel].get()->getConnectionHandlerTrialSpace(),ieti,m_parSpaceHandler[spaceLevel],settings_.precondIterations,settings_.precondTol,false,false,maxEigenvalue,outputIETI_iter);
                                    maxEigenvalue =dynamic_cast<IETIAdapterMPI<T>*>(precond.get())->getMaxEigenvalue();
                                }
                                    break;
                                case Base::MG:
                                {
                                    size_t nLevels = spaceLevel+1;
                                    if(nLevels>1)
                                    {
                                        std::vector<typename gsDistributedOperator<T>::Ptr> linOpsD (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOps (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsR (nLevels-1);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsP (nLevels-1);
                                        for(size_t l=0; l<nLevels;++l)
                                        {
                                            linOpsD[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                            linOps[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                        }
                                        for(size_t l =0; l<nLevels-1;++l)
                                        {
                                            linOpsR[l] = m_restriction[l];
                                            linOpsP[l] = m_prolongation[l];
                                        }

                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);

                                        typename gsMultiGridOp<T>::Ptr mg = gsMultiGridOp<T>::make(linOps,linOpsP,linOpsR,cs );
                                        mg->setOptions(settings_.MGOptions);

                                        for (int l = mg->numLevels() == 1 ? 0 : 1; l < mg->numLevels(); ++l)
                                        {
                                            switch( settings_.smoother ) {
                                            case Smoother_::Richardson:                            mg->setSmoother(l, makeRichardsonOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::Jacobi:                                mg->setSmoother(l, makeJacobiOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::GaussSeidel:                           mg->setSmoother(l, makeGaussSeidelOp(mg->matrix(l))); break;
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD:
                                            {
                                                // std::vector<typename gsLinearOperator<T>::Ptr> localSmoothers = makeSubspaceCorrectedMassSmootherOperatorsDirichlet(settings_.m_spatialAssemblers[i]->multiBasis(),settings_.damping);
                                                GISMO_NO_IMPLEMENTATION;
                                                /*
                                             mg->setSmoother(
                                                         i,
                                                         gsParallelAdditivePreconditionerOp<T>::make(
                                                             linOpsD[i],
                                                             setupPiecewisePreconditioner<T>(
                                                                 *linOpsD[i], //NeedsParallelOperator!!!
                                                                 *m_parSpaceHandler[i],
                                                                 give(localSmoothers),
                                                                 settings_.m_spatialAssemblers[i]->patches(),
                                                                 settings_.m_spatialAssemblers[i]->multiBasis(),
                                                                 settings_.m_spatialAssemblers[i]->pde().bc(),
                                                                 settings_.m_spatialAssemblers[i]->options(),
                                                                 myPatches_,
                                                                 comm_
                                                                 ),
                                                             settings_.outerDamping
                                                             )
                                                         );
                                                         */
                                                break;
                                            }
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:
                                            {
                                                index_t degree = settings_.m_spatialAssemblers[l]->multiBasis().maxCwiseDegree();
                                                const index_t dim = settings_.m_spatialAssemblers[l]->multiBasis().dim();
                                                gsMultiBasis<T> proc_loc_mb;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                    proc_loc_mb.addBasis(settings_.m_spatialAssemblers[l]->multiBasis()[myPatches_[k]].clone().release());

                                                const index_t proc_loc_dofs = (*m_parSpaceHandler[l]).localSize();

                                                gsDofMapper dm; // Global mapper
                                                settings_.m_spatialAssemblers[l]->multiBasis().getMapper(
                                                            (dirichlet::strategy)settings_.m_spatialAssemblers[l]->options().askInt("DirichletStrategy",11),
                                                            (iFace    ::strategy)settings_.m_spatialAssemblers[l]->options().askInt("InterfaceStrategy", 1),
                                                            settings_.m_spatialAssemblers[l]->pde().bc(),
                                                            dm,
                                                            0
                                                            );

                                                // Global => Processor local
                                                gsMatrix<index_t> glob2loc;
                                                {
                                                    gsMatrix<index_t> inp(proc_loc_dofs,1);
                                                    inp.col(0) = gsVector<index_t>::LinSpaced(proc_loc_dofs,0,proc_loc_dofs-1);
                                                    glob2loc.setZero((*m_parSpaceHandler[l]).globalSize(),1);
                                                    (*m_parSpaceHandler[l]).addLocalVectorToGlobal( inp, glob2loc );
                                                }

                                                std::vector< gsVector<index_t> > proc_loc_maps;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                {
                                                    const index_t sz = proc_loc_mb[k].size();
                                                    index_t kk = myPatches_[k];
                                                    gsVector<index_t> local(sz);
                                                    for (index_t j=0; j<sz; ++j)
                                                    {
                                                        local[j] = dm.is_free(j,kk) ? glob2loc(dm.index(j, kk),0) : -1u;
                                                    }
                                                    proc_loc_maps.push_back(give(local));
                                                }
                                                std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<> > > > pieces =
                                                        constructPieces( proc_loc_mb, proc_loc_maps, proc_loc_dofs );

                                                std::vector< typename gsLinearOperator<T>::Ptr > localSmoothers;
                                                std::vector< gsSparseMatrix<T> > smootherTransfers;
                                                const index_t nrPieces = pieces.size();

                                                real_t h = 1;
                                                for ( size_t jj=0; jj<settings_.m_spatialAssemblers[l]->multiBasis().nBases(); ++jj)
                                                {
                                                    const gsGeometry<T>& patch = settings_.m_spatialAssemblers[l]->patches().patch(jj);
                                                    real_t diam =(patch.coefAtCorner(boxCorner::getFirst(dim)) - patch.coefAtCorner(boxCorner::getLast(dim))).norm();
                                                    h = std::min(h,settings_.m_spatialAssemblers[l]->multiBasis()[jj].getMinCellLength()*diam);
                                                }
                                                for ( index_t dd = 0; dd<nrPieces; ++dd )
                                                {
                                                    gsBoundaryConditions<T> localbc;
                                                    for( index_t ps=0; ps < 2*dim; ++ps )
                                                        localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

                                                    const index_t sz = pieces[dd].size();
                                                    for ( index_t jj=0; jj<sz; ++jj)
                                                    {
                                                        if (dd == dim)
                                                        {
                                                            localSmoothers.push_back( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc,alpha));
                                                        }
                                                        else if (dd > 1)
                                                        {
                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            localSmoothers.push_back( gsScaledOp<T>::make( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc, (scalingM+(alpha*math::pow(2,dim-dd)*math::pow(h/(2*degree+1),dim-dd) )/scalingK)), 1/scalingK ) );
                                                        }
                                                        else if (dd == 1)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );

                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            gsSparseMatrix<T> M, K;
                                                            assembleParameterMass(*pieces[dd][jj].first, M);
                                                            assembleParameterStiffness(*pieces[dd][jj].first, K);
                                                            M = M.block( 1, 1, M.rows()-2, M.cols()-2 );
                                                            K = K.block( 1, 1, K.rows()-2, K.cols()-2 );
                                                            gsSparseMatrix<T> mat = scalingK * K + (scalingM+(alpha*math::pow(2,dim-dd)*math::pow( h/(2*degree+1), dim-dd ))) * M;

                                                            localSmoothers.push_back( makeSparseCholeskySolver(mat) );
                                                        }
                                                        else if (dd == 0)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );
                                                            const index_t sz = pieces[dd][jj].second.rows();
                                                            const real_t scalingFactor = math::pow(2,dim) * dim * math::pow( h/(1+2*degree), dim-1 ) * degree*degree/(h*(2*degree-1))+ math::pow(h/(2*degree+1),dim)*math::pow(2,dim)*alpha;
                                                            localSmoothers.push_back( gsScaledOp<T>::make( gsIdentityOp<T>::make(sz), 1/scalingFactor) );
                                                        }

                                                        smootherTransfers.push_back( give( pieces[dd][jj].second ) );
                                                    }
                                                }

                                                mg->setSmoother(
                                                            l,
                                                            gsParallelAdditivePreconditionerOp<T>::make(
                                                                linOpsD[l],
                                                                give(smootherTransfers),
                                                                give(localSmoothers),
                                                                settings_.outerDamping
                                                                )
                                                            );
                                                break;
                                            }


                                            }
                                        }

                                        precond = gsParallelPreconditionerOp<T>::make(mg);
                                      //  gsInfo<<"precond size: "<<precond->rows()<<"  -  "<<"mg size: "<<mg->rows()<<"\n";
                                    }
                                    else
                                    {
                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);
                                        precond = gsParallelOperator<T>::make(m_Kx[spaceLevel]->getConnectionHandlerTrialSpace(),m_Kx[spaceLevel]->getConnectionHandlerTestSpace(),cs);
                                      //  gsInfo<<"precond size: "<<precond->rows()<<"  -  "<<"cs size: "<<cs->rows()<<"\n";
                                    }

                                    break;
                                }
                                };

                                APrec[i]=ComplexSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],(*UT)(i,i),precond,comm_,settings_.approxIterates,settings_.tol,true);
                                AEx[i]=ComplexSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],(*UT)(i,i),precond,comm_,100,tolExactInnerSolv,false);


                                /*
                                typename gsMatrix<T>::EigenSolver a;
                                gsMatrix<T> locMat,locMat_, Prec,col,col_;
                                gsMatrix<T> K,M;
                                gsMatrix<T> MK;
                                if(slice == 0)
                                {
                                    gsMatrix<T> id = gsMatrix<T>::Identity(m_parSpaceHandler[spaceLevel]->globalSize(),m_parSpaceHandler[spaceLevel]->globalSize());
                                    m_parSpaceHandler[spaceLevel]->extractLocalVector(id,locMat_);
                                    precond->distribute(locMat_,locMat);
                                    for(int c=0; c<locMat.cols();++c)
                                    {
                                        precond->apply(locMat.col(c),col);
                                        precond->distribute(col,col_);
                                        locMat.col(c) = col_;
                                    }
                                    m_parSpaceHandler[spaceLevel]->buildGlobalVector(locMat,Prec);
                                    a.compute(Prec);
                                    T symmNorm = (Prec-Prec.transpose()).norm();
                                    if(comm_.rank()==0) gsInfo<<"\t\t MG: eigenvalues for lambda: "<<(*UT)(i,i)<<" -> "<<a.eigenvalues().real().minCoeff()<<" - "<<a.eigenvalues().real().maxCoeff()<<" || symm: "<<symmNorm<<" \n";


                                    K= m_Kx[spaceLevel]->toDenseMatrix(*m_parSpaceHandler[spaceLevel]);
                                    M= m_Mx[spaceLevel]->toDenseMatrix(*m_parSpaceHandler[spaceLevel]);

                                    MK = K + ((*UT)(i,i)).real()*M;
                                    symmNorm = (MK-MK.transpose()).norm();
                                    a.compute(MK);
                                    if(comm_.rank()==0) gsInfo<<"\t\t KM eigenvalues for lambda: "<<(*UT)(i,i)<<" -> "<<a.eigenvalues().real().minCoeff()<<" - "<<a.eigenvalues().real().maxCoeff()<<" || symm: "<<symmNorm<<" \n";
                                    //a.compute(Prec*MK);
                                    //if(comm_.rank()==0) gsInfo<<"\t\t MG-KM: eigenvalues for lambda: "<<(*UT)(i,i)<<" -> "<<a.eigenvalues().minCoeff()<<" - "<<a.eigenvalues().maxCoeff()<<"\n";
                                }
                                */

                                for(int j=i+1;j<nt;++j)
                                {
                                    if(std::abs((*UT)(i,i) - std::conj((*UT)(j,j)))<2.e-10)
                                    {
                                        APrec[j]=ComplexSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],(*UT)(j,j),precond,comm_,settings_.approxIterates,settings_.tol,true);
                                        AEx[j]=ComplexSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],(*UT)(j,j),precond,comm_,100,tolExactInnerSolv,false);

                                        /*
                                        if(false && slice == 0)
                                        {
                                            MK = K + ((*UT)(j,j)).real()*M;
                                            a.compute(MK);
                                            if(comm_.rank()==0) gsInfo<<"\t\t KM eigenvalues for lambda: "<<(*UT)(j,j)<<" -> "<<a.eigenvalues().real().minCoeff()<<" - "<<a.eigenvalues().real().maxCoeff()<<"\n";
                                            //a.compute(Prec*MK);
                                            //if(comm_.rank()==0) gsInfo<<"\t\t MG-KM: eigenvalues for lambda: "<<(*UT)(i,i)<<" -> "<<a.eigenvalues().minCoeff()<<" - "<<a.eigenvalues().maxCoeff()<<"\n";
                                        }
                                        */
                                        break;
                                    }
                                }
                            }
                        }

                        precA[level][slice] = gsDeComplexify<T>::make(  gsProductOp<std::complex<T> >::make(
                                                                            gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(Pinv),gsIdentityOp<std::complex<T> >::make(nx)),
                                                                            Base::BlockTriangularSolver::make(APrec,m_Mx[spaceLevel],UT),
                                                                            gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(Q),gsIdentityOp<std::complex<T> >::make(nx))));

                        LuOfA[level][slice] = gsDeComplexify<T>::make( gsProductOp<std::complex<T> >::make(
                                                                           gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(Pinv),gsIdentityOp<std::complex<T> >::make(nx)),
                                                                           Base::BlockTriangularSolver::make(AEx,m_Mx[spaceLevel],UT),
                                                                           gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(Q),gsIdentityOp<std::complex<T> >::make(nx))));

                        /*
                        gsMatrix<T> K,M,Ainv;
                        K = m_Kx[spaceLevel]->toDenseMatrix(*m_parSpaceHandler[spaceLevel]);
                        M = m_Mx[spaceLevel]->toDenseMatrix(*m_parSpaceHandler[spaceLevel]);

                        for(int i=0; i<nt;++i)
                        {
                            gsMatrix<std::complex<T> > AEXinvI;
                            AEx[i]->toMatrix(AEXinvI);
                            gsInfo<<"Ai^{-1}: \n"<<AEXinvI<<"\n\n";
                            gsMatrix<std::complex<T> > Ai = K + (*UT)(i,i)* M;
                            gsInfo<<"Ai^{-1}*Ai on "<<i<<"  : "<<(Ai*AEXinvI - gsMatrix<std::complex<T> >::Identity(Ai.rows(),Ai.rows())).array().abs().maxCoeff()<<"\n";
                        }

                        gsMatrix<T> A = gsMatrix<T>(*data_[level].Kt[slice]).kron(M) + gsMatrix<T>(*data_[level].Mt[slice]).kron(K);
                        gsInfo<<"A=\n"<<A<<"\n\n";
                        typename Base::BlockTriangularSolver::Ptr solver =  Base::BlockTriangularSolver::make(AEx,m_Mx[spaceLevel],UT);
                        gsMatrix<std::complex<T> > result;
                        solver->toMatrix(result);
                        gsInfo<<"UT=\n"<<*UT<<"\n\n";
                        gsInfo<<"Tri=\n"<<result<<"\n\n";
                        LuOfA[level][slice]->toMatrix(Ainv);
                        gsInfo<<"Ainv=\n"<<Ainv<<"\n\n";
                        gsInfo<<"A^{-1}*A  : "<<(A*Ainv - gsMatrix<std::complex<T> >::Identity(A.rows(),A.rows())).array().abs().maxCoeff()<<"\n";

                        abort();
                        */
                        //     typename BlockTriangularSolver::Ptr solver =  BlockTriangularSolver::make(APrec,data_[level].Mx[slice],UT);
                        //       gsMatrix<std::complex<T> > result;
                        //      solver->toMatrix(result);

                        //        gsInfo<<"A: \n"<< A2.toDense()<<"\n\n";
                        //         gsInfo<<"A^{-1}: \n"<< A2.toDense().inverse()<<"\n\n";
                        //         gsInfo<<"pA: \n"<< result.inverse()<<"\n\n";
                        //          gsInfo<<"pA^{-1}: \n"<< result<<"\n\n";
                    }
                        break;
                    case Base::REALSCHUR:
                    {
                        gsMatrix<T> Minv = Mt.inverse();
                        Eigen::RealSchur<Eigen::Matrix<T,Dynamic,Dynamic> > solv(Minv*Kt, true);
                        const gsMatrix<T >& TT = solv.matrixT();

                        gsMatrix<T> givens = gsMatrix<T>::Identity(nt,nt);

                        int nBlocks = 0;
                        for(int i=0;i<nt;++i)
                        {
                            if((i==nt-1 && math::abs(TT(i,i-1)) <complexBound) || (i!=nt-1  && math::abs(TT(i+1,i)) <complexBound) )
                                nBlocks++;
                            else
                            {
                                T t = (TT(i,i+1) + TT(i+1,i))/(TT(i,i)-TT(i+1,i+1));
                                T t_s = math::min(-t - math::sqrt(t*t + 1),-t + math::sqrt(t*t + 1));
                                T c = T(1)/math::sqrt(1+t_s*t_s);
                                T s = t_s*c;
                                givens(i,i) = c;
                                givens(i+1,i+1) = c;
                                givens(i,i+1)= s;
                                givens(i+1,i) = -s;
                                nBlocks++;
                                i++;
                            }

                        }

                        typename gsMatrix<T >::Ptr Q = (gsMatrix<T>(solv.matrixU()*givens)).moveToPtr();
                        typename gsMatrix<T >::Ptr Pinv =  (gsMatrix<T>(Q->transpose()*Minv)).moveToPtr();

                        typename gsMatrix<T >::Ptr UT = (gsMatrix<T>(givens.transpose()*TT*givens)).moveToPtr();

                        std::vector<typename gsLinearOperator<T >::Ptr> APrec(nBlocks);
                        std::vector<typename gsLinearOperator<T >::Ptr> AEx(nBlocks);
                        int idx = 0;
                        for(int i=0;i<nBlocks;++i)
                        {
                            GISMO_ASSERT((*UT)(idx, idx)>=0, "The real part of the Eigenvalues becomes negative, choose a smaller theta");
                            //(K + (alpha+|beta|)M)^{-1}
                            if((i==nBlocks-1 && idx +1 ==nt) || (i!=nBlocks-1  && math::abs((*UT)(idx+1,idx)) <complexBound) )
                            {
                                T alpha = (*UT)(idx, idx);
                                switch(settings_.spaceSolver){
                                case Base::DIRECT:
                                    GISMO_ERROR("Direct Solvers not supported in space-parallel mode");
                                    //APrec[i]=gsSolverOp<sparseLLTfact>::make(((*UT)(idx,idx)*(*data_[level].Mx[slice]) + *data_[level].Kx[slice]));
                                    //AEx[i]=APrec[i];
                                    break;
                                case Base::IETI:
                                {
                                    typename gsIETIAssemblerMPI<T>::Ptr ieti = memory::make_shared(new gsIETIAssemblerMPI<T>(*settings_.m_spatialAssemblers[spaceLevel],myPatches_,comm_));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.reserve(data_[level].PatchMx[slice].size());
                                    for(size_t npi = 0; npi<myPatches_.size();++npi)
                                    {
                                        size_t np = myPatches_[npi];
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    alpha * (*data_[level].PatchMx[slice][np])
                                                    + *data_[level].PatchKx[slice][np]);
                                        mats.push_back(give(mat));
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();

                                    APrec[i] = IETIAdapterMPI<T>::make(m_Kx[spaceLevel].get()->getConnectionHandlerTrialSpace(),ieti,m_parSpaceHandler[spaceLevel],settings_.precondIterations,settings_.precondTol,false,false,maxEigenvalue,outputIETI_iter);
                                    AEx[i] =IETIAdapterMPI<T>::make(m_Kx[spaceLevel].get()->getConnectionHandlerTrialSpace(),ieti,m_parSpaceHandler[spaceLevel],100,tolExactInnerSolv,false,false,-1,outputIETI_iter);
                                    maxEigenvalue  =dynamic_cast<IETIAdapterMPI<T>*>(APrec[i].get())->getMaxEigenvalue();

                                }
                                case Base::MG:
                                {
                                    size_t nLevels = spaceLevel+1;
                                    if(nLevels>1)
                                    {
                                        std::vector<typename gsDistributedOperator<T>::Ptr> linOpsD (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOps (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsR (nLevels-1);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsP (nLevels-1);
                                        for(size_t l=0; l<nLevels;++l)
                                        {
                                            linOpsD[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                            linOps[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                        }
                                        for(size_t l =0; l<nLevels-1;++l)
                                        {
                                            linOpsR[l] = m_restriction[l];
                                            linOpsP[l] = m_prolongation[l];
                                        }

                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);

                                        typename gsMultiGridOp<T>::Ptr mg = gsMultiGridOp<T>::make(linOps,linOpsP,linOpsR,cs);
                                        mg->setOptions(settings_.MGOptions);

                                        for (int l = mg->numLevels() == 1 ? 0 : 1; l < mg->numLevels(); ++l)
                                        {
                                            switch( settings_.smoother ) {
                                            case Smoother_::Richardson:                            mg->setSmoother(l, makeRichardsonOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::Jacobi:                                mg->setSmoother(l, makeJacobiOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::GaussSeidel:                           mg->setSmoother(l, makeGaussSeidelOp(mg->matrix(l))); break;
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD:
                                            {
                                                GISMO_NO_IMPLEMENTATION;
                                                //  std::vector<typename gsLinearOperator<T>::Ptr> localSmoothers = makeSubspaceCorrectedMassSmootherOperatorsDirichlet(settings_.m_spatialAssemblers[i]->multiBasis(),settings_.damping);
                                                /*
                                             mg->setSmoother(
                                                         i,
                                                         gsParallelAdditivePreconditionerOp<T>::make(
                                                             linOpsD[i],
                                                             setupPiecewisePreconditioner<T>(
                                                                 *linOpsD[i], //NeedsParallelOperator!!!
                                                                 *m_parSpaceHandler[i],
                                                                 give(localSmoothers),
                                                                 settings_.m_spatialAssemblers[i]->patches(),
                                                                 settings_.m_spatialAssemblers[i]->multiBasis(),
                                                                 settings_.m_spatialAssemblers[i]->pde().bc(),
                                                                 settings_.m_spatialAssemblers[i]->options(),
                                                                 myPatches_,
                                                                 comm_
                                                                 ),
                                                             settings_.outerDamping
                                                             )
                                                         );
                                                         */
                                                break;
                                            }
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:
                                            {
                                                index_t degree = settings_.m_spatialAssemblers[l]->multiBasis().maxCwiseDegree();
                                                const index_t dim = settings_.m_spatialAssemblers[l]->multiBasis().dim();
                                                gsMultiBasis<T> proc_loc_mb;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                    proc_loc_mb.addBasis(settings_.m_spatialAssemblers[l]->multiBasis()[myPatches_[k]].clone().release());

                                                const index_t proc_loc_dofs = (*m_parSpaceHandler[l]).localSize();

                                                gsDofMapper dm; // Global mapper
                                                settings_.m_spatialAssemblers[l]->multiBasis().getMapper(
                                                            (dirichlet::strategy)settings_.m_spatialAssemblers[l]->options().askInt("DirichletStrategy",11),
                                                            (iFace    ::strategy)settings_.m_spatialAssemblers[l]->options().askInt("InterfaceStrategy", 1),
                                                            settings_.m_spatialAssemblers[l]->pde().bc(),
                                                            dm,
                                                            0
                                                            );

                                                // Global => Processor local
                                                gsMatrix<index_t> glob2loc;
                                                {
                                                    gsMatrix<index_t> inp(proc_loc_dofs,1);
                                                    inp.col(0) = gsVector<index_t>::LinSpaced(proc_loc_dofs,0,proc_loc_dofs-1);
                                                    glob2loc.setZero((*m_parSpaceHandler[l]).globalSize(),1);
                                                    (*m_parSpaceHandler[l]).addLocalVectorToGlobal( inp, glob2loc );
                                                }

                                                std::vector< gsVector<index_t> > proc_loc_maps;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                {
                                                    const index_t sz = proc_loc_mb[k].size();
                                                    index_t kk = myPatches_[k];
                                                    gsVector<index_t> local(sz);
                                                    for (index_t j=0; j<sz; ++j)
                                                    {
                                                        local[j] = dm.is_free(j,kk) ? glob2loc(dm.index(j, kk),0) : -1u;
                                                    }
                                                    proc_loc_maps.push_back(give(local));
                                                }
                                                std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<> > > > pieces =
                                                        constructPieces( proc_loc_mb, proc_loc_maps, proc_loc_dofs );

                                                std::vector< typename gsLinearOperator<T>::Ptr > localSmoothers;
                                                std::vector< gsSparseMatrix<T> > smootherTransfers;
                                                const index_t nrPieces = pieces.size();

                                                real_t h = 1;
                                                for ( size_t jj=0; jj<settings_.m_spatialAssemblers[l]->multiBasis().nBases(); ++jj)
                                                {
                                                    const gsGeometry<T>& patch = settings_.m_spatialAssemblers[l]->patches().patch(jj);
                                                    real_t diam =(patch.coefAtCorner(boxCorner::getFirst(dim)) - patch.coefAtCorner(boxCorner::getLast(dim))).norm();
                                                    h = std::min(h,settings_.m_spatialAssemblers[l]->multiBasis()[jj].getMinCellLength()*diam);
                                                }
                                                for ( index_t dd = 0; dd<nrPieces; ++dd )
                                                {
                                                    gsBoundaryConditions<T> localbc;
                                                    for( index_t ps=0; ps < 2*dim; ++ps )
                                                        localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

                                                    const index_t sz = pieces[dd].size();
                                                    for ( index_t jj=0; jj<sz; ++jj)
                                                    {
                                                        if (dd == dim)
                                                        {
                                                            localSmoothers.push_back( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc,alpha) );
                                                        }
                                                        else if (dd > 1)
                                                        {
                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            localSmoothers.push_back( gsScaledOp<T>::make( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc, (scalingM+alpha*math::pow(2,dim-dd)*math::pow(h/(2*degree+1),dim-dd) )/scalingK), 1/scalingK ) );
                                                        }
                                                        else if (dd == 1)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );

                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            gsSparseMatrix<T> M, K;
                                                            assembleParameterMass(*pieces[dd][jj].first, M);
                                                            assembleParameterStiffness(*pieces[dd][jj].first, K);
                                                            M = M.block( 1, 1, M.rows()-2, M.cols()-2 );
                                                            K = K.block( 1, 1, K.rows()-2, K.cols()-2 );
                                                            gsSparseMatrix<T> mat = scalingK * K + scalingM * M + alpha*scalingK/(dim-dd) *M;

                                                            localSmoothers.push_back( makeSparseCholeskySolver(mat) );
                                                        }
                                                        else if (dd == 0)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );
                                                            const index_t sz = pieces[dd][jj].second.rows();
                                                            const real_t scalingFactor = std::pow(2,dim) * dim * std::pow( h/(1+2*degree), dim-1 ) * degree*degree/(h*(2*degree-1))+ alpha*std::pow(2,dim)*std::pow( h/(1+2*degree), dim);
                                                            localSmoothers.push_back( gsScaledOp<T>::make( gsIdentityOp<T>::make(sz), 1/scalingFactor) );
                                                        }

                                                        smootherTransfers.push_back( give( pieces[dd][jj].second ) );
                                                    }
                                                }

                                                mg->setSmoother(
                                                            l,
                                                            gsParallelAdditivePreconditionerOp<T>::make(
                                                                linOpsD[l],
                                                                give(smootherTransfers),
                                                                give(localSmoothers),
                                                                settings_.outerDamping
                                                                )
                                                            );
                                                break;
                                            }
                                            }
                                        }

                                        typename gsParallelPreconditionerOp<T>::Ptr mg_prec = gsParallelPreconditionerOp<T>::make(mg);

                                        APrec[i]  = mg_prec;
                                        gsOptionList options;
                                        options.addReal("Iterations","number of maximal Iterations", 100);
                                        options.addReal("Tolerance", "tolerance", tolExactInnerSolv);
                                        AEx[i] = gsDistributedIterativeSolverOp<gsParallelCG<T> >::make(m_Kx.back()->getConnectionHandlerTestSpace(),linOpsD.back(), mg_prec,options);
                                      //  gsInfo<<"AEx[i] size: "<<AEx[i]->rows()<<"mg-size: "<<mg_prec<<"\n";
                                    }
                                    else
                                    {
                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);
                                       // gsInfo<<"cs size: "<<cs->rows()<<"\n";
                                        APrec[i]  =  cs;
                                        AEx[i] = cs;

                                    }
                                }
                                    break;
                                default:
                                    GISMO_NO_IMPLEMENTATION;
                                };
                            }
                            else
                            {
                                T alpha = (*UT)(idx, idx) + math::sqrt(math::abs((*UT)(idx + 1, idx) * (*UT)(idx, idx + 1)));
                                //the abs is needed here, since the sign is not changed at this point
                                typename gsDistributedOperator<T>::Ptr precond;
                                switch(settings_.spaceSolver){
                                case Base::DIRECT:
                                    GISMO_ERROR("Direct Solvers not supported in space-parallel mode");
                                    // precond= gsSolverOp<sparseLLTfact>::make((gsSparseMatrix<T>(*data_[level].Kx[slice]+alpha**data_[level].Mx[slice])));
                                    break;
                                case Base::IETI:
                                {
                                    typename gsIETIAssemblerMPI<T>::Ptr ieti = memory::make_shared(new gsIETIAssemblerMPI<T>(*settings_.m_spatialAssemblers[spaceLevel],myPatches_,comm_));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.reserve(data_[level].PatchMx[slice].size());
                                    for(size_t npi = 0; npi<myPatches_.size();++npi)
                                    {
                                        size_t np = myPatches_[npi];
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    alpha * (*data_[level].PatchMx[slice][np]) + *data_[level].PatchKx[slice][np]);
                                        mats.push_back(give(mat));
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();

                                    precond = IETIAdapterMPI<T>::make(m_Kx[spaceLevel].get()->getConnectionHandlerTrialSpace(),ieti,m_parSpaceHandler[spaceLevel],settings_.precondIterations,settings_.precondTol,false,false,maxEigenvalue,outputIETI_iter);
                                    maxEigenvalue  =dynamic_cast<IETIAdapterMPI<T>*>(precond.get())->getMaxEigenvalue();
                                }
                                    break;
                                case Base::MG:
                                {

                                    size_t nLevels = spaceLevel+1;
                                    if(nLevels>1)
                                    {
                                        std::vector<typename gsDistributedOperator<T>::Ptr> linOpsD (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOps (nLevels);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsR (nLevels-1);
                                        std::vector<typename gsLinearOperator<T>::Ptr> linOpsP (nLevels-1);
                                        for(size_t l=0; l<nLevels;++l)
                                        {
                                            linOpsD[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                            linOps[l] = gsDistributedSumOp<T>::make(m_Kx[l],gsDistributedScaledOp<T>::make(m_Mx[l],alpha));
                                        }
                                        for(size_t l =0; l<nLevels-1;++l)
                                        {
                                            linOpsR[l] = m_restriction[l];
                                            linOpsP[l] = m_prolongation[l];
                                        }

                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);

                                        typename gsMultiGridOp<T>::Ptr mg = gsMultiGridOp<T>::make(linOps,linOpsP,linOpsR,cs );
                                        mg->setOptions(settings_.MGOptions);

                                        for (int l = mg->numLevels() == 1 ? 0 : 1; l < mg->numLevels(); ++l)
                                        {
                                            switch( settings_.smoother ) {
                                            case Smoother_::Richardson:                            mg->setSmoother(l, makeRichardsonOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::Jacobi:                                mg->setSmoother(l, makeJacobiOp(mg->matrix(l),settings_.damping*settings_.outerDamping)); break;
                                            case Smoother_::GaussSeidel:                           mg->setSmoother(l, makeGaussSeidelOp(mg->matrix(l))); break;
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD:
                                            {
                                                // std::vector<typename gsLinearOperator<T>::Ptr> localSmoothers = makeSubspaceCorrectedMassSmootherOperatorsDirichlet(settings_.m_spatialAssemblers[i]->multiBasis(),settings_.damping);
                                                GISMO_NO_IMPLEMENTATION;
                                                /*
                                             mg->setSmoother(
                                                         i,
                                                         gsParallelAdditivePreconditionerOp<T>::make(
                                                             linOpsD[i],
                                                             setupPiecewisePreconditioner<T>(
                                                                 *linOpsD[i], //NeedsParallelOperator!!!
                                                                 *m_parSpaceHandler[i],
                                                                 give(localSmoothers),
                                                                 settings_.m_spatialAssemblers[i]->patches(),
                                                                 settings_.m_spatialAssemblers[i]->multiBasis(),
                                                                 settings_.m_spatialAssemblers[i]->pde().bc(),
                                                                 settings_.m_spatialAssemblers[i]->options(),
                                                                 myPatches_,
                                                                 comm_
                                                                 ),
                                                             settings_.outerDamping
                                                             )
                                                         );
                                                         */
                                                break;
                                            }
                                            case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:
                                            {
                                                index_t degree = settings_.m_spatialAssemblers[l]->multiBasis().maxCwiseDegree();
                                                const index_t dim = settings_.m_spatialAssemblers[l]->multiBasis().dim();
                                                gsMultiBasis<T> proc_loc_mb;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                    proc_loc_mb.addBasis(settings_.m_spatialAssemblers[l]->multiBasis()[myPatches_[k]].clone().release());

                                                const index_t proc_loc_dofs = (*m_parSpaceHandler[l]).localSize();

                                                gsDofMapper dm; // Global mapper
                                                settings_.m_spatialAssemblers[l]->multiBasis().getMapper(
                                                            (dirichlet::strategy)settings_.m_spatialAssemblers[l]->options().askInt("DirichletStrategy",11),
                                                            (iFace    ::strategy)settings_.m_spatialAssemblers[l]->options().askInt("InterfaceStrategy", 1),
                                                            settings_.m_spatialAssemblers[l]->pde().bc(),
                                                            dm,
                                                            0
                                                            );

                                                // Global => Processor local
                                                gsMatrix<index_t> glob2loc;
                                                {
                                                    gsMatrix<index_t> inp(proc_loc_dofs,1);
                                                    inp.col(0) = gsVector<index_t>::LinSpaced(proc_loc_dofs,0,proc_loc_dofs-1);
                                                    glob2loc.setZero((*m_parSpaceHandler[l]).globalSize(),1);
                                                    (*m_parSpaceHandler[l]).addLocalVectorToGlobal( inp, glob2loc );
                                                }

                                                std::vector< gsVector<index_t> > proc_loc_maps;
                                                for (size_t k=0;k<myPatches_.size();++k)
                                                {
                                                    const index_t sz = proc_loc_mb[k].size();
                                                    index_t kk = myPatches_[k];
                                                    gsVector<index_t> local(sz);
                                                    for (index_t j=0; j<sz; ++j)
                                                    {
                                                        local[j] = dm.is_free(j,kk) ? glob2loc(dm.index(j, kk),0) : -1u;
                                                    }
                                                    proc_loc_maps.push_back(give(local));
                                                }
                                                std::vector< std::vector< std::pair< typename gsBasis<T>::Ptr, gsSparseMatrix<> > > > pieces =
                                                        constructPieces( proc_loc_mb, proc_loc_maps, proc_loc_dofs );

                                                std::vector< typename gsLinearOperator<T>::Ptr > localSmoothers;
                                                std::vector< gsSparseMatrix<T> > smootherTransfers;
                                                const index_t nrPieces = pieces.size();

                                                real_t h = 1;
                                                for ( size_t jj=0; jj<settings_.m_spatialAssemblers[l]->multiBasis().nBases(); ++jj)
                                                {
                                                    const gsGeometry<T>& patch = settings_.m_spatialAssemblers[l]->patches().patch(jj);
                                                    real_t diam =(patch.coefAtCorner(boxCorner::getFirst(dim)) - patch.coefAtCorner(boxCorner::getLast(dim))).norm();
                                                    h = std::min(h,settings_.m_spatialAssemblers[l]->multiBasis()[jj].getMinCellLength()*diam);
                                                }

                                                for ( index_t dd = 0; dd<nrPieces; ++dd )
                                                {
                                                    gsBoundaryConditions<T> localbc;
                                                    for( index_t ps=0; ps < 2*dim; ++ps )
                                                        localbc.addCondition( 0, 1+ps, condition_type::dirichlet, NULL );

                                                    const index_t sz = pieces[dd].size();
                                                    for ( index_t jj=0; jj<sz; ++jj)
                                                    {
                                                        if (dd == dim)
                                                        {
                                                            localSmoothers.push_back( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc,alpha ));
                                                        }
                                                        else if (dd > 1)
                                                        {
                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            localSmoothers.push_back( gsScaledOp<T>::make( makeSubspaceCorrectedMassSmootherOperator(*pieces[dd][jj].first, settings_.damping, localbc, (scalingM+(alpha*math::pow(2,dim-dd)*math::pow(h/(2*degree+1),dim-dd) )/scalingK)), 1/scalingK ) );
                                                        }
                                                        else if (dd == 1)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );

                                                            const real_t scalingK = math::pow(2,dim-dd) * (dim-dd) * math::pow( h/(2*degree+1), dim-dd );
                                                            const real_t scalingM = math::pow(2,dim-dd) * math::pow( h/(2*degree+1), dim-dd-1 ) * degree*degree/(h*(2*degree-1));
                                                            gsSparseMatrix<T> M, K;
                                                            assembleParameterMass(*pieces[dd][jj].first, M);
                                                            assembleParameterStiffness(*pieces[dd][jj].first, K);
                                                            M = M.block( 1, 1, M.rows()-2, M.cols()-2 );
                                                            K = K.block( 1, 1, K.rows()-2, K.cols()-2 );
                                                            gsSparseMatrix<T> mat = scalingK * K + scalingM * M + alpha*scalingK/(dim-dd) *M;

                                                            localSmoothers.push_back( makeSparseCholeskySolver(mat) );
                                                        }
                                                        else if (dd == 0)
                                                        {
                                                            //gsSparseMatrix<> localMat = pieces[dd][j].second * mg->matrix(i) * pieces[dd][j].second.transpose();
                                                            //localSmoothers.push_back( makeSparseCholeskySolver(localMat) );
                                                            const index_t sz = pieces[dd][jj].second.rows();
                                                            const real_t scalingFactor = std::pow(2,dim) * dim * std::pow( h/(1+2*degree), dim-1 ) * degree*degree/(h*(2*degree-1))+ alpha*std::pow(2,dim)*std::pow( h/(1+2*degree), dim);
                                                            localSmoothers.push_back( gsScaledOp<T>::make( gsIdentityOp<T>::make(sz), 1/scalingFactor) );
                                                        }

                                                        smootherTransfers.push_back( give( pieces[dd][jj].second ) );
                                                    }
                                                }

                                                mg->setSmoother(
                                                            l,
                                                            gsParallelAdditivePreconditionerOp<T>::make(
                                                                linOpsD[l],
                                                                give(smootherTransfers),
                                                                give(localSmoothers),
                                                                settings_.outerDamping
                                                                )
                                                            );
                                                break;
                                            }

                                            }
                                        }

                                        precond = gsParallelPreconditionerOp<T>::make(mg);
                                      //  gsInfo<<"precond size: "<<precond->rows()<<"  -  "<<"mg size: "<<mg->rows()<<"\n";
                                    }
                                    else
                                    {
                                        std::vector<typename  gsMatrixOp<gsSparseMatrix<T> >::Ptr> c_ops(data_[level].PatchMx[slice].size());
                                        for(size_t k=0; k<myPatches_.size();++k)
                                        {
                                            size_t np = myPatches_[k];
                                            c_ops[np]=makeMatrixOp(gsSparseMatrix<T>(m_coarseKx[np]->matrix()+m_coarseMx[np]->matrix()*alpha).moveToPtr());
                                        }

                                        typename gsCoarseSolverAdapter<T>::Ptr cs = constructCoarseSolver<T>(c_ops,m_parSpaceHandler.front(),myPatches_,m_coarseLocMappers,m_coarseGlobMapper);
                                        precond = gsParallelOperator<T>::make(m_Kx[spaceLevel]->getConnectionHandlerTrialSpace(),m_Kx[spaceLevel]->getConnectionHandlerTestSpace(),cs);
                                      //  gsInfo<<"precond size: "<<precond->rows()<<"  -  "<<"cs size: "<<cs->rows()<<"\n";
                                    }
                                    break;
                                }
                                default:
                                    GISMO_NO_IMPLEMENTATION;
                                };

                                APrec[i]=RealSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],(*UT)(idx,idx),(*UT)(idx+1,idx),(*UT)(idx,idx+1),precond,settings_.approxIterates,settings_.tol);
                                AEx[i]=RealSolverMPI::make(m_Kx[spaceLevel],m_Mx[spaceLevel],(*UT)(idx,idx),(*UT)(idx+1,idx),(*UT)(idx,idx+1),precond,100,tolExactInnerSolv);
                               // gsInfo<<"RealSolverMPI size: "<<AEx[i]->rows() <<" x "<<AEx[i]->cols()<<"\n";

                                idx++;
                            }
                            idx++;
                        }

                        /*{
                            typename Base::QuasiBlockTriangularSolver::uPtr test = Base::QuasiBlockTriangularSolver::make(APrec,m_Mx[spaceLevel],UT);
                            gsInfo<<"QuasiBlockTriangularSolver size: "<<test->rows() <<" x "<<test->cols()<<"\n";
                        }*/
                        precA[level][slice] =gsProductOp<T>::make(   gsKroneckerOp<T>::make(makeMatrixOp(Pinv),gsIdentityOp<T>::make(nx)),
                                                                     Base::QuasiBlockTriangularSolver::make(APrec,m_Mx[spaceLevel],UT),
                                                                     gsKroneckerOp<T>::make(makeMatrixOp(Q),gsIdentityOp<T>::make(nx)));

                        LuOfA[level][slice]= gsProductOp<T>::make(  gsKroneckerOp<T>::make(makeMatrixOp(Pinv),gsIdentityOp<T>::make(nx)),
                                                                    Base::QuasiBlockTriangularSolver::make(AEx,m_Mx[spaceLevel],UT),
                                                                    gsKroneckerOp<T>::make(makeMatrixOp(Q),gsIdentityOp<T>::make(nx)));
                        break;
                    }
                    default:
                        gsInfo<<"No other decompostition implemented.";
                        GISMO_NO_IMPLEMENTATION;
                        break;
                    };
                }
            }
            //   if(output)std::cout<<"finished calc of Ainv on level "<<level<< std::endl;+
            break;
        }
            //  */
        default:
            gsInfo<<"no other solver implemented, choose ("<<Base::FULLSPACETIME<<") for direct solver or ("<<Base::SPATIALWISE<<") for direct solver utilizing the tensor product structure.";
            GISMO_NO_IMPLEMENTATION;
            break;

        };

    }
}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::setMGCoarseMatrices(std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr > coarseKx,
                                                      std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr > coarseMx,
                                                      std::vector<gsDofMapper> coarseLocMappers,
                                                      gsDofMapper coarseGlobMapper)
{
    m_coarseKx = coarseKx;
    m_coarseMx = coarseMx;
    m_coarseLocMappers = coarseLocMappers;
    m_coarseGlobMapper = coarseGlobMapper;
}


template<typename T>
void gsSpaceParallelHSTSlapOp<T>::extractAndReorderRhs(gsMatrix<T>& rhs) const
{
    gsMatrix<T> tempRhs;

    m_parHandler->extractLocalVector(rhs,tempRhs,comm_.rank());
    //gsInfo<<"rank: "<<m_parHandler->getComm().rank()<<" local rhs: "<<tempRhs.rows()<<" and local size handler: "<<m_parHandler->localSize()<<"\n";
    rhs = tempRhs;
}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::buildSliceWiseSolution(std::vector<gsMatrix<T> >& sol, const std::vector<int> & sliceIdx)const
{
    gsMatrix<T> glob;
    for(size_t ii=0; ii<sliceIdx.size();++ii)
    {
        glob.setZero(m_parHandler->globalSize(),sol[ii].cols());

        dynamic_cast<gsParallelOperator<T>*>(A.back()[ii].get())->getConnectionHandlerTrialSpace()->distributeAccumulatedVector(sol[ii]);
        m_parHandler->buildGlobalVector(sol[ii],glob);
        sol[ii] = glob;
        //   gsInfo<<"glob: \n"<<glob.transpose()<<"\n\n"<<"sol: \n"<<sol[ii].transpose()<<"\n\n\n";
    }
}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::buildSliceWiseRHS(std::vector<gsMatrix<T> >& rhs, const std::vector<int> & sliceIdx)const
{
    gsMatrix<T> glob;
    for(size_t ii=0; ii<sliceIdx.size();++ii)
    {
        glob.setZero(m_parHandler->globalSize(),rhs[ii].cols());
        m_parHandler->buildGlobalVector(rhs[ii],glob);

        rhs[ii] = glob;
    }
}


template<typename T>
void gsSpaceParallelHSTSlapOp<T>::mult(const Vector &u, Vector &f, int timeLevel, int timeStep, real_t sign) const
{
    gsMatrix<T> result;
    A[timeLevels_[timeLevel]][timeStep]->apply(u, result);
    f +=sign*result;

    /*
    gsMatrix<T> globVec;
  //  m_parHandler.buildGlobalVector(u,globVec);
  //  if(gsMpi::worldRank()==0) gsInfo<<"before apply A: \n"<<globVec.transpose()<<"\n\n";
    m_parHandler.buildGlobalVector(result,globVec);
    if(gsMpi::worldRank()==0) gsInfo<<"after apply A: \n"<<globVec.transpose()<<"\n\n";
    m_parHandler.buildGlobalVector(f,globVec);
    if(gsMpi::worldRank()==0) gsInfo<<"after apply A +f: \n"<<globVec.transpose()<<"\n\n";
//*/
}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::calcInitialVector(const Vector &u, Vector &f, int timeLevel, int timeStep) const
{
    //GISMO_ASSERT(timeStep>0, "timeStep is 0, no B is defined here");
    Vector res;
    B[timeLevels_[timeLevel]][timeStep]->apply(u,res);
    f+=res;

    /*
    gsMatrix<T> globVec;
    m_parHandler.buildGlobalVector(f,globVec);
    if(gsMpi::worldRank()==0) gsInfo<<"after apply B: \n"<<globVec.transpose()<<"\n\n";
//*/
}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::TimeRestriction1(const Vector &u_fine, Vector &u_coarse, int coarseTimeLevel, bool setZero) const
{
    GISMO_UNUSED(coarseTimeLevel);
    /*
    gsMatrix<T> globVec;
    m_parHandler.buildGlobalVector(u_fine,globVec);
    if(gsMpi::worldRank()==0) gsInfo<<"before restr1: \n"<<globVec.transpose()<<"\n\n";
//*/
    gsKroneckerOp<T> kron(makeMatrixOp(timeRestr1_),gsIdentityOp<T>::make((index_t)u_fine.rows()/timeRestr1_.cols()));

    if(setZero) u_coarse.setZero();
    kron.apply(u_fine,m_temp);
    u_coarse += m_temp;

    /*
    m_parHandler.buildGlobalVector(m_temp,globVec);
    if(gsMpi::worldRank()==0) gsInfo<<"after restr1: \n"<<globVec.transpose()<<"\n\n";
//*/
}

//TODO: identicalSlice
template<typename T>
void gsSpaceParallelHSTSlapOp<T>::TimeRestriction2(const Vector &u_fine, Vector &u_coarse, int coarseTimeLevel, bool setZero) const
{
    GISMO_UNUSED(coarseTimeLevel);
    /*
    gsMatrix<T> globVec;
    m_parHandler.buildGlobalVector(u_fine,globVec);
    if(gsMpi::worldRank()==0) gsInfo<<"before restr2: \n"<<globVec.transpose()<<"\n\n";
//*/
    gsKroneckerOp<T> kron(makeMatrixOp(timeRestr2_),gsIdentityOp<T>::make((index_t)u_fine.rows()/timeRestr2_.cols()));

    if(setZero) u_coarse.setZero();
    kron.apply(u_fine,m_temp);
    u_coarse += m_temp;

    /*
    m_parHandler.buildGlobalVector(m_temp,globVec);
    if(gsMpi::worldRank()==0) gsInfo<<"after restr2: \n"<<globVec.transpose()<<"\n\n";
    //*/
}

//TODO: identicalSlice
template<typename T>
void gsSpaceParallelHSTSlapOp<T>::TimeProlongation1(const Vector &u_coarse, Vector &u_fine, int coarseTimeLevel, bool setZero) const
{
    GISMO_UNUSED(coarseTimeLevel);
    //gsKroneckerOp<T> kron(makeMatrixOp(timeRestr1_.transpose()),gsIdentityOp<T>::make((index_t)m_parHandler.localSize()/timeRestr1_.cols()));
    gsKroneckerOp<T> kron(makeMatrixOp(timeRestr1_.transpose()),gsIdentityOp<T>::make((index_t)u_coarse.rows()/timeRestr1_.rows()));
    if(setZero) u_fine.setZero();
    kron.apply(u_coarse,m_temp);
    u_fine += m_temp;

    /*
    gsMatrix<T> globVec;
    gsMatrix<T> uFine = u_fine;
    dynamic_cast<gsParallelOperator<T>*>(A[coarseTimeLevel].front().get())->getConnectionHandlerTestSpace()->distributeAccumulatedVector(uFine);
    m_parHandler.buildGlobalVector(uFine,globVec);
    if(gsMpi::worldRank()==0) gsInfo<<"after prol1: \n"<<globVec.transpose()<<"\n\n";
    //*/
}

//TODO:identicalSlice
template<typename T>
void gsSpaceParallelHSTSlapOp<T>::TimeProlongation2(const Vector &u_coarse, Vector &u_fine, int coarseTimeLevel, bool setZero) const
{
    GISMO_UNUSED(coarseTimeLevel);
    //gsKroneckerOp<T> kron(makeMatrixOp(timeRestr2_.transpose()),gsIdentityOp<T>::make((index_t)m_parHandler.localSize()/timeRestr2_.cols()));
    gsKroneckerOp<T> kron(makeMatrixOp(timeRestr2_.transpose()),gsIdentityOp<T>::make((index_t)u_coarse.rows()/timeRestr2_.rows()));
    if(setZero) u_fine.setZero();
    kron.apply(u_coarse,m_temp);
    u_fine += m_temp;
    /*
    gsMatrix<T> globVec;
    gsMatrix<T> uFine = u_fine;
    dynamic_cast<gsParallelOperator<T>*>(A[coarseTimeLevel].front().get())->getConnectionHandlerTestSpace()->distributeAccumulatedVector(uFine);
    m_parHandler.buildGlobalVector(uFine,globVec);
    if(gsMpi::worldRank()==0) gsInfo<<"after prol2: \n"<<globVec.transpose()<<"\n\n";
    //*/
}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::accumulate(const Vector & v,Vector & v_acc, MPI_Comm spaceComm) const
{
    GISMO_UNUSED(spaceComm);
    v_acc = v;
    dynamic_cast<gsParallelOperator<T>* >(A.back().front().get())->getConnectionHandlerTestSpace()->accumulateDistributedVector(v_acc);
}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::distribute(const Vector & v,Vector & v_dist, MPI_Comm spaceComm) const
{
    GISMO_UNUSED(spaceComm);
    v_dist = v;
    dynamic_cast<gsParallelOperator<T>* >(A.back().front().get())->getConnectionHandlerTestSpace()->distributeAccumulatedVector(v_dist);
}



template<typename T>
gsSpaceParallelHSTSlapOp<T>::ComplexSolverMPI::ComplexSolverMPI(typename gsParallelOperator<T>::Ptr K, typename gsParallelOperator<T>::Ptr M, const std::complex<T>& lambda, typename gsDistributedOperator<T>::Ptr precond, gsMpiComm comm, int max_iter, real_t tol, bool prec, std::string outputFilename) : m_K(K), m_M(M), m_lambda(lambda)
{
    m_filename = outputFilename;
    m_max_iter = max_iter;
    m_tol = tol;

    m_inp.resize(2*cols(),1);
    m_out.resize(2*rows(),1);

    blockOP = gsDistributedBlockOp<T>::make(2,2);
    blockOP->addOperator(0,0, gsDistributedSumOp<T>::make(K,gsDistributedScaledOp<T>::make(M,lambda.real())));
    blockOP->addOperator(1,1, gsDistributedSumOp<T>::make(gsDistributedScaledOp<T>::make(K,-1),gsDistributedScaledOp<T>::make(M,-lambda.real())));
    blockOP->addOperator(0,1, gsDistributedScaledOp<T>::make(M,lambda.imag()));
    blockOP->addOperator(1,0, gsDistributedScaledOp<T>::make(M,lambda.imag()));

    blockPrec = gsDistributedBlockOp<T>::make(2,2);
    blockPrec->addOperator(0,0, precond);
    blockPrec->addOperator(1,1, precond);
    itSolver = gsParallelMinRes<T>::make(blockOP,blockPrec,comm);
    itSolver->setMaxIterations(max_iter);
    itSolver->setTolerance(tol);
    m_print= true;
    //if(prec)
    itSolver->setInexactResidual(true);
}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::ComplexSolverMPI::apply(const gsMatrix<std::complex<T> > & input, gsMatrix<std::complex<T> > & x) const
{
    x.resize(rows(),input.cols());
    m_out.setZero(2*rows(),input.cols());
    m_inp.resize(2*rows(),input.cols());
    m_inp.topRows(rows()) = input.real();
    m_inp.bottomRows(rows()) = input.imag();

    if(input.cols()>1)
    {
        gsMatrix<T > temp = gsMatrix<T>::Zero(m_inp.rows(),1);
        for(int i=0; i< input.cols();++i)
        {
            itSolver->solve(m_inp.col(i), temp );
            m_out.col(i) = temp;
        }
    }
    else
        itSolver->solve(m_inp, m_out);

    x.real() = m_out.topRows(rows());
    x.imag() = -m_out.bottomRows(rows());
    if(m_print)
    {
        std::fstream out(m_filename.c_str(),std::ofstream::out | std::ofstream::app);
        out<<"Complex solver numIT: "<<itSolver->iterations()<<"  - res: "<<itSolver->error()<<"  - size ("<<input.rows()<<")\n";
        gsInfo<<"Complex solver numIT: "<<itSolver->iterations()<<"  - res: "<<itSolver->error()<<" - size ("<<input.rows()<<")\n";
    }
}


template<typename T>
gsSpaceParallelHSTSlapOp<T>::RealSolverMPI::RealSolverMPI(typename gsParallelOperator<T>::Ptr K, typename gsParallelOperator<T>::Ptr M, T a, T b, T c , typename gsDistributedOperator<T>::Ptr precond, int max_iter, real_t tol, std::string outputFilename) : m_K(K), m_M(M), m_a(a), m_b(b), m_c(c)
{
    // m_MinRes = MinRes;
    m_max_iter = max_iter;
    m_tol = tol;
    m_filename = outputFilename;

    m_inp.resize(cols(),1);
    m_out.resize(rows(),1);

    blockOP = gsDistributedBlockOp<T>::make(2,2);
    /*
    blockOP->addOperator(0,0, gsSumOp<T>::make(makeMatrixOp(K),gsScaledOp<T>::make(makeMatrixOp(M),m_a)));
    blockOP->addOperator(1,1, gsSumOp<T>::make(gsScaledOp<T>::make(makeMatrixOp(K),-1),gsScaledOp<T>::make(makeMatrixOp(M),-m_a)));
    blockOP->addOperator(0,1, gsScaledOp<T>::make(makeMatrixOp(M),-m_c));
    blockOP->addOperator(1,0, gsScaledOp<T>::make(makeMatrixOp(M),m_b));

    blockPrec = gsBlockOp<T>::make(2,2);
    blockPrec->addOperator(0,0, precond);
    blockPrec->addOperator(1,1, precond);
*/

    blockOP->addOperator(0,0, gsDistributedSumOp<T>::make(gsDistributedScaledOp<T>::make(K,math::abs(m_b)),gsDistributedScaledOp<T>::make(M,math::abs(m_b)*m_a)));
    blockOP->addOperator(1,1, gsDistributedSumOp<T>::make(gsDistributedScaledOp<T>::make(K,-math::abs(m_c)),gsDistributedScaledOp<T>::make(M,-math::abs(m_c)*m_a)));
    blockOP->addOperator(0,1, gsDistributedScaledOp<T>::make(M,-math::abs(m_b)*m_c));
    blockOP->addOperator(1,0, gsDistributedScaledOp<T>::make(M,math::abs(m_c)*m_b));

    blockPrec = gsDistributedBlockOp<T>::make(2,2);
    blockPrec->addOperator(0,0,  gsDistributedScaledOp<T>::make(precond, T(1)/math::abs(m_b)));
    blockPrec->addOperator(1,1, gsDistributedScaledOp<T>::make(precond,T(1)/math::abs(m_c)));

    /*
    gsMatrix<T> big, pre;
    blockOP->toMatrix(big);
    blockPrec->toMatrix(pre);
    gsInfo<<"mat: \n"<<big<<std::endl;
    gsInfo<<"prec: \n"<<pre<<std::endl;
*/ //blockPrec
    //   itSolver = memory::make_shared(new gsGMRes<>(blockOP,blockPrec));
    itSolver = memory::make_shared(new gsParallelMinRes<T>(blockOP,blockPrec,K->getConnectionHandlerTestSpace()->getComm()));

    itSolver->setMaxIterations(max_iter);
    itSolver->setTolerance(tol);
    itSolver->setInexactResidual(true);
    m_print= true;

}

template<typename T>
void gsSpaceParallelHSTSlapOp<T>::RealSolverMPI::apply(const gsMatrix<T > & input, gsMatrix<T> & x) const
{
    x.resize(rows(),input.cols());
    x.setZero();
    gsMatrix<T>& inp = const_cast<gsMatrix<T>&>(input);
    inp.topRows(blockOP->getOperator(0,0)->rows())*=math::abs(m_b);
    inp.bottomRows(blockOP->getOperator(1,1)->rows())*=math::abs(m_c);
    if(inp.cols()>1)
        for(int i=0; i< inp.cols();++i)
        {
            itSolver->solve(inp.col(i), m_temp );
            x.col(i) = m_temp;
        }
    else
        itSolver->solve(inp, x);
    //Revert changes in input
    inp.topRows(blockOP->getOperator(0,0)->rows())/=math::abs(m_b);
    inp.bottomRows(blockOP->getOperator(1,1)->rows())/=math::abs(m_c);

    //  gsMatrix<T> result;
    //  blockOP->apply(x,result);
    //T error = (input-result).norm()/input.norm();
    x.bottomRows(blockOP->getOperator(1,1)->rows())*=-1;
    if(m_print)
    {
        std::fstream out(m_filename.c_str(),std::ofstream::out | std::ofstream::app);
        out<<"Real solver numIT: "<<itSolver->iterations()<<"  - res: "<<itSolver->error()<<"\n";
        gsInfo<<"\tReal solver numIT: "<<itSolver->iterations()<<"  - res: "<<itSolver->error()<<" - size ("<<input.rows()<<")\n";
    }
}


} // namespace gismo
#endif

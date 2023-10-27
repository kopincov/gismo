/**  gsMGSpaceTimeAdapter.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    Created on:  2017-11-23
*/


#include <gsAssembler/gsMGSpaceTimeAdapter.h>
#include <gsSolver/gsSumOp.h>
#include <gsMultiGrid/gsMultiGrid.h>
#include <gsIETI/gsIETIAdapter.h>
#include <gsSolver/gsComplexify.h>
#include <fstream>

namespace gismo {


template<typename T>
gsHeatSTSlapOperator<T>::HeatSTSlapSettings::HeatSTSlapSettings(gsOptionList opt)
{
    gsOptionList dOpt = defaultOptions();
    dOpt.update(opt);
    sliceSolver = dOpt.getInt("Solver");
    spaceSolver = dOpt.getInt("DecompositionSolver");
    timeDecomp = dOpt.getInt("DecompositionMethod");
    approxIterates = dOpt.getInt("Iterations");
    tol = dOpt.getReal("Tolerance");
    precondIterations = dOpt.getInt("PreconditionIterations");
    precondTol = dOpt.getReal("PreconditionTol");
    smoother = Smoother_::type(dOpt.getInt("Smoother"));
    damping = dOpt.getReal("damping");
    outerDamping = dOpt.getReal("outerDamping");
    IETIOptions = dOpt.getGroup("IETI");
    MGOptions = dOpt.getGroup("MG");
    //  IETIOptions.setSwitch("NoAssemblyOfRHS",true);
    if(spaceSolver==IETI && sliceSolver == FULLSPACETIME)
        IETIOptions.setSwitch("NonSymmetric",true);

    if (damping < 0)
    {
        // Here, one could add appropriate parameter-choice rules (depending on the degree)
        switch (smoother)
        {
        case Smoother_::Richardson:                                     damping = 0.80; break;
        case Smoother_::Jacobi:                                         damping = 0.80; break;
        case Smoother_::GaussSeidel:                                    break;
        case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD:       damping = 0.09; break;
        case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:      damping = 0.09; break;
        }
    }
}

template<typename T>
gsHeatSTSlapOperator<T>::HeatSTSlapSettings::HeatSTSlapSettings()
{
    gsOptionList dOpt = defaultOptions();
    sliceSolver = dOpt.getInt("Solver");
    spaceSolver = dOpt.getInt("DecompositionSolver");
    timeDecomp = dOpt.getInt("DecompositionMethod");
    approxIterates = dOpt.getInt("Iterations");
    tol = dOpt.getReal("Tolerance");
    precondIterations = dOpt.getInt("PreconditionIterations");
    precondTol = dOpt.getReal("PreconditionTol");
    smoother = Smoother_::type(dOpt.getInt("Smoother"));
    damping = dOpt.getReal("damping");
    outerDamping = dOpt.getReal("outerDamping");
    IETIOptions = dOpt.getGroup("IETI");
    MGOptions = dOpt.getGroup("MG");

    //  IETIOptions.setSwitch("NoAssemblyOfRHS",true);
    if(spaceSolver==IETI && sliceSolver == FULLSPACETIME)
        IETIOptions.setSwitch("NonSymmetric",true);

    if (damping < 0)
    {
        // Here, one could add appropriate parameter-choice rules (depending on the degree)
        switch (smoother)
        {
        case Smoother_::Richardson:                                     damping = 0.80; break;
        case Smoother_::Jacobi:                                         damping = 0.80; break;
        case Smoother_::GaussSeidel:                                    break;
        case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD:       damping = 0.09; break;
        case Smoother_::MassRichardsonSubspaceCorrectionAdditiveMPDDD2:      damping = 0.09; break;
        }
    }

}

template<typename T>
gsOptionList gsHeatSTSlapOperator<T>::HeatSTSlapSettings::defaultOptions()
{
    gsOptionList opt;
    opt = gsIETIAssembler<T>::defaultOptions().wrapIntoGroup("IETI");
    opt.update(gsMultiGridOp<T>::defaultOptions().wrapIntoGroup("MG"),gsOptionList::addIfUnknown);
    opt.addReal("damping",            "Damping factor for the smoother (handed over to smoother)",-1);
    opt.addReal("outerDamping",       "Damping factor for the smoother (globally)",1.0);
    opt.addInt("Smoother","chosen smoother: (0) Richardson, (1) Jacobi, (2) GS, (3) SubCorr_MPDDD, (4) SubCorr_MPDDD2 (only this works in parallel)", 4);
    opt.addInt("Solver", "chosen solver for A block, (0) FullSpaceTime, (1) spatialwise slolver,  (2) MG (not implemented) (3) GMRES (not implemented)",0);
    opt.addInt("DecompositionSolver", "chosen method for spatial problems, (0) direct solver, (1) IETI", 0);
    opt.addInt("DecompositionMethod", "chosed matrix decomposition for A, (0) diagonalization (may not work), (1) Complex Schur, (2) Real Schur",2);
    opt.addInt("Iterations","max number of iterations for the approximative solve",10);
    opt.addReal("Tolerance", "tolerance for the iterative solver for the subproblems in A",1.e-4);
    opt.addInt("PreconditionIterations", "max number of iteraions for preconditioner in approximative solve",5);
    opt.addReal("PreconditionTol", "tolerance for the preconditioner in approximative solve", 1.e-4);
    return opt;
}


template<typename T>
gsHeatSTSlapOperator<T>::gsHeatSTSlapOperator(real_t h, real_t tau, int NTimeLevels)  { h_ = h; tau_ = tau; NTimeLevels_ = NTimeLevels; ndofs_.resize(1);}

//Matrices are already in Tensor form
template<typename T>
gsHeatSTSlapOperator<T>::gsHeatSTSlapOperator(const std::vector<int> &TimeLevels, int NTimeLevels,
                                              int spaceLevels,
                                              int spaceLevel,
                                              std::vector<std::vector<int> > myTimeSlices,
                                              std::vector<STData > data,
                                              std::vector<typename gsSparseMatrix<T,RowMajor>::Ptr >& transfer,
                                              gsMatrix<T> timeRestr1,
                                              gsMatrix<T> timeRestr2,
                                              real_t tau,
                                              real_t h,
                                              typename gsHeatSTSlapOperator<T>::HeatSTSlapSettings settings, bool output)
    : settings_(settings), myTimeSlices(myTimeSlices), timeRestr1_(timeRestr1), timeRestr2_(timeRestr2), TimeLevelsInp_(TimeLevels), transferInp_(transfer), output_(output)
{
    ndofs_.resize(1,0);
    this->spaceLevel = spaceLevel;
    if(TimeLevelsInp_.size()==0) return;

    //    real_t t00 = omp_get_wtime();
    if(output) std::cout << "initialize block multigrid solver..." << std::endl;

    // Mh_ = Mh.back();
    // m_ = Mh.back()->m(); n_ = Mh.back()->n();

    spaceLevels_ = spaceLevels;//Mh.size();

    NTimeLevels_ = NTimeLevels;
    //TODO: identicalSlice
    h_ = h;  //mesh.getMaxH(spaceLevels_-1);

    //compute the time matrices
    //TODO: identicalSlice
    tau_ = tau;
    //  order_ = order;
    //  MatTimeK_.resize(TimeLevels.size());
    //   MatTimeM_.resize(TimeLevels.size());
    //   MatTimeT0_.resize(TimeLevels.size());
    timeLevels_.resize(NTimeLevels_); for(int i=0; i<NTimeLevels_; i++) timeLevels_[i] = -1;
    for(unsigned int i=0; i<TimeLevelsInp_.size(); i++)
    {
        timeLevels_[TimeLevelsInp_[i]] = i;
        // int level = TimeLevels[i];
        //   real_t tauL = real_t(pow(2,NTimeLevels_-level-1))*tau_;
        //   assembleParabolicMatricies(tauL, order_, MatTimeK_[i], MatTimeM_[i], MatTimeT0_[i]);
    }
    //    if(output)std::cout << "did timelevels..." << std::endl;

    // if(TimeLevels.size()>0) blocksize_ = MatTimeM_[0].m();
    //  else blocksize_ = order+1;

    //compute the time restriction matrix
    /*        Identity id;
        int dimension = 1;



        Vector p0(1); p0[0] = 0.0;
        Vector p1(1); p1[0] = 0.5;
        Vector p2(1); p2[0] = 1.0;

        std::vector<Vector> pELCoarse(0); pELCoarse.push_back(p0); pELCoarse.push_back(p2);
        std::vector<Vector> pELFine1(0); pELFine1.push_back(p0); pELFine1.push_back(p1);
        std::vector<Vector> pELFine2(0); pELFine2.push_back(p1); pELFine2.push_back(p2);

        ElementMapping elCoarse(pELCoarse, 1, true);
        ElementMapping elFine1(pELFine1, 1, true);
        ElementMapping elFine2(pELFine2, 1, true);

        //a temporary lookUpTable
        std::vector<std::vector<std::vector<int> > > lookUpTable(0);
        CGFE fel(dimension, order_, lookUpTable);

        Matrix mass1;
        Matrix mass2;
        calcElmatMass(elFine1, fel, id, mass1);
        calcElmatMass(elFine2, fel, id, mass2);

        calcElmatMass(elFine1, elCoarse, elFine1, fel, fel, id, timeRestr1_); timeRestr1_ = timeRestr1_ | inverse(mass1);
        calcElmatMass(elFine2, elCoarse, elFine2, fel, fel, id, timeRestr2_); timeRestr2_ = timeRestr2_ | inverse(mass2);
*/
    this->data_ = data;
}

//Matrices are already in Tensor form
template<typename T>
gsHeatSTSlapOperator<T>::gsHeatSTSlapOperator(const std::vector<int> &TimeLevels, int NTimeLevels,
                                              int spaceLevels,
                                              int spaceLevel,
                                              std::vector<std::vector<int> > myTimeSlices,
                                              std::vector<STData > data,
                                              gsMatrix<T> timeRestr1,
                                              gsMatrix<T> timeRestr2,
                                              real_t tau,
                                              real_t h,
                                              typename gsHeatSTSlapOperator<T>::HeatSTSlapSettings settings, bool output)
    : settings_(settings), myTimeSlices(myTimeSlices), timeRestr1_(timeRestr1), timeRestr2_(timeRestr2), TimeLevelsInp_(TimeLevels), output_(output)
{
    ndofs_.resize(1,0);
    this->spaceLevel = spaceLevel;
    if(TimeLevelsInp_.size()==0) return;

    if(output) std::cout << "initialize block multigrid solver..." << std::endl;
    spaceLevels_ = spaceLevels;

    NTimeLevels_ = NTimeLevels;
    h_ = h;
    tau_ = tau;
    timeLevels_.resize(NTimeLevels_); for(int i=0; i<NTimeLevels_; i++) timeLevels_[i] = -1;
    for(unsigned int i=0; i<TimeLevelsInp_.size(); i++)
    {
        timeLevels_[TimeLevelsInp_[i]] = i;
    }
    this->data_ = data;
}

template<typename T>
void gsHeatSTSlapOperator<T>::init()
{
    A.resize(TimeLevelsInp_.size());
    B.resize(TimeLevelsInp_.size());
    bool isEmpty = true;
    for(unsigned int level=0; level<TimeLevelsInp_.size(); level++)
    {
        size_t size = myTimeSlices[TimeLevelsInp_[level]].size();

        A[level].resize(size);
        B[level].resize(size);
        if(size>0) isEmpty = false;
        for(size_t slice = 0; slice<size;++slice) //data[level].Kt.size()
        {
            int tSlice = myTimeSlices[TimeLevelsInp_[level]][slice];
            if(output_)std::cout<<" , level "<<level<<", slice "<<slice <<"  calculating A,B, mapped slice "<<tSlice << std::endl;
            A[level][slice] = gsSumOp<T>::make(
                        gsKroneckerOp<T>::make(makeMatrixOp(data_[level].Kt[tSlice]),makeMatrixOp(data_[level].Mx[slice])),
                        gsKroneckerOp<T>::make(makeMatrixOp(data_[level].Mt[tSlice]),makeMatrixOp(data_[level].Kx[slice])));
            if(myTimeSlices[TimeLevelsInp_[level]][slice]>0)
                B[level][slice] = gsKroneckerOp<T>::make(makeMatrixOp(data_[level].Nt[tSlice]),makeMatrixOp(data_[level].MWx[slice]));
            //Changed slice to tSlice
        }
    }
    if(output_)std::cout<<"finished calc of A,B "<< std::endl;
    if(!isEmpty)ndofs_.resize(data_.back().Mt.size());
    for(size_t sl = 0; sl< ndofs_.size();++sl)
        ndofs_[sl] = timeRestr1_.cols()* (spaceLevel == 0? transferInp_.front()->cols(): transferInp_[spaceLevel-1]->rows());//A.back()[sl]->rows(); //THIS MUST WORK IN ALL CASES!
    if(output_)std::cout<<"finished ndofs "<< std::endl;
    if(ndofs_.front() == 0) return;


    /*
        A_.resize(TimeLevels.size());
        int blocks = blocksize_;
        localMatrix coeff(blocks, blocks);
        localMatrix filling(blocks, blocks);

        for(unsigned int level=0; level<TimeLevels.size(); level++)
        {
            A_[level].resize(spaceLevels_);

            //Stiffness-localMatrix
            coeff = 0.0; filling = 0.0;
            for(int k=0; k<blocksize_; k++) for(int l=0; l<blocksize_; l++)
            {
                coeff(k,l) = MatTimeM_[level](k,l);
                filling(k,l) = 1.0;
            }
            for(int xLevel=0; xLevel<spaceLevels_; xLevel++) A_[level][xLevel].addMatrix(Kh[xLevel], coeff, filling, false);


            //Mass-localMatrix
            coeff = 0.0; filling = 0.0;
            for(int k=0; k<blocksize_; k++) for(int l=0; l<blocksize_; l++)
            {
                coeff(k,l) = MatTimeK_[level](k,l);
                filling(k,l) = 1.0;
            }
            for(int xLevel=0; xLevel<spaceLevels_; xLevel++) A_[level][xLevel].addMatrix(Mh[xLevel], coeff, filling, false);


        }
*/


    //ndofs_ = A_.back().back().n();

    int nRestr = transferInp_.size();
    spaceRestr_.resize(nRestr);
    for(int i=0; i<nRestr; i++) { spaceRestr_[i] = transferInp_[i]; }
    // if(output)std::cout<<"finished restr "<< std::endl;
    transferSpace.resize(nRestr);
    transferTSpace.resize(nRestr);
    /*
        for(unsigned int step=0; step<spaceRestr_.size(); step++)
        {
            coeff = 0.0; filling = 0.0;
            for(int k=0; k<blocksize_; k++)
            {
                coeff(k,k) = 1.0;
                filling(k,k) = 1.0;
            }
            spaceRestr_[step].addMatrix(restr_[step], coeff, filling, false);
        }
*/
    for(unsigned int step=0; step<transferSpace.size(); step++)
    {
        //spaceRestr_[step] is actually a prolongation :D
        transferTSpace[step] = gsKroneckerOp<T>::make(gsIdentityOp<T>::make(timeRestr1_.cols()),makeMatrixOp(spaceRestr_[step]));
        transferSpace[step] = gsKroneckerOp<T>::make(gsIdentityOp<T>::make(timeRestr1_.cols()),makeMatrixOp(spaceRestr_[step]->transpose()));
    }
    //   if(output)std::cout<<"finished Space transfer "<< std::endl;
    real_t maxEigenvalue = -1; //for IETI
    switch(settings_.sliceSolver)
    {
    case FULLSPACETIME:
    {
        precA.resize(TimeLevelsInp_.size());
        LuOfA.resize(TimeLevelsInp_.size());
        for(unsigned int level=0; level<TimeLevelsInp_.size(); level++)
        {
            precA[level].resize(myTimeSlices[TimeLevelsInp_[level]].size());
            LuOfA[level].resize(myTimeSlices[TimeLevelsInp_[level]].size());
            for(size_t slice = 0; slice<myTimeSlices[TimeLevelsInp_[level]].size();++slice)
            {

                switch(settings_.spaceSolver){
                case DIRECT:
                {
                    gsSparseMatrix<T> Mat, temp;
                    temp=data_[level].Kt[myTimeSlices[TimeLevelsInp_[level]][slice]]->kron(*data_[level].Mx[slice]);
                    Mat=data_[level].Mt[myTimeSlices[TimeLevelsInp_[level]][slice]]->kron(*data_[level].Kx[slice]);
                    Mat+=temp;

                    precA[level][slice] = gsSolverOp<sparseLUfact>::make(Mat);
                    LuOfA[level][slice]  = precA[level][slice];
                    break;
                }
                case IETI:
                {
                    typename gsIETIAssembler<T>::Ptr ieti = memory::make_shared(new gsIETIAssembler<T>(*settings_.m_STAssemblers[spaceLevel]));

                    ieti->setOptions(settings_.IETIOptions);
                    ieti->init();
                    std::vector<gsSparseMatrix<T> > mats;
                    mats.reserve(data_[level].PatchMx[slice].size());
                    gsSparseMatrix<T> Mat, temp;
                    for(size_t np=0; np<data_[level].PatchMx[slice].size();++np)
                    {
                        temp=data_[level].Kt[myTimeSlices[TimeLevelsInp_[level]][slice]]->kron(*data_[level].PatchMx[slice][np]);
                        Mat=data_[level].Mt[myTimeSlices[TimeLevelsInp_[level]][slice]]->kron(*data_[level].PatchKx[slice][np]);
                        Mat+=temp;
                        mats.push_back(give(Mat));
                    }
                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                    ieti->assemble();

                    LuOfA[level][slice] = IETIAdapter<T>::make(ieti,100,1.e-10,false,true);
                    precA[level][slice]  = IETIAdapter<T>::make(ieti,settings_.precondIterations,settings_.precondTol,true,true,maxEigenvalue);
                    maxEigenvalue = dynamic_cast<IETIAdapter<T>*>(precA[level][slice].get())->getMaxEigenvalue();
                    break;
                }
                }

                //   if(output)std::cout<<"finished calc of Ainv on level "<<level<< std::endl;
            }

        }
        break;
    case SPATIALWISE:
        {
            precA.resize(TimeLevelsInp_.size());
            LuOfA.resize(TimeLevelsInp_.size());
            for(unsigned int level=0; level<TimeLevelsInp_.size(); level++)
            {
                precA[level].resize(myTimeSlices[TimeLevelsInp_[level]].size());
                LuOfA[level].resize(myTimeSlices[TimeLevelsInp_[level]].size());
                for(size_t slice = 0; slice<myTimeSlices[TimeLevelsInp_[level]].size();++slice)
                {

                    int nx = data_[level].Mx[slice]->cols();
                    int nt = data_[level].Mt[myTimeSlices[TimeLevelsInp_[level]][slice]]->cols();

                    gsMatrix<T> Kt = *data_[level].Kt[myTimeSlices[TimeLevelsInp_[level]][slice]];
                    gsMatrix<T> Mt = *data_[level].Mt[myTimeSlices[TimeLevelsInp_[level]][slice]];

                    switch(settings_.timeDecomp)
                    {
                    case DIAGONALIZATION:
                    {
                        Eigen::GeneralizedEigenSolver<Eigen::Matrix<T,Dynamic,Dynamic> > solv;
                        solv.compute(Kt,Mt,true);
                        typename Eigen::GeneralizedEigenSolver<Eigen::Matrix<T,Dynamic,Dynamic> >::EigenvalueType D = solv.eigenvalues();
                        typename gsMatrix<std::complex<T> >::Ptr X = gsMatrix<std::complex<T> > (solv.eigenvectors()).moveToPtr();
                        typename gsMatrix<std::complex<T> >::Ptr VinvTr = gsMatrix<std::complex<T> > ((Mt*(*X))).moveToPtr();
                        //  gsInfo<<"Kt: \n"<<Kt<<"\n\n"<<"Mt: \n"<<Mt<<"\n\n";
                        //   gsInfo<<"X: \n"<<solv.eigenvectors()<<"\n\n"<<"D: \n"<<D<<"\n\n"<<"V^-T: \n"<<*VinvTr<<"\n\n";

                        typename gsBlockOp<std::complex<T> >::Ptr blockPrec = gsBlockOp<std::complex<T> >::make(nt,nt);
                        typename gsBlockOp<std::complex<T> >::Ptr blockEx = gsBlockOp<std::complex<T> >::make(nt,nt);
                        for(int i=0;i<nt;++i)
                        {
                            //(K + (alpha+|beta|)M)^{-1}

                            if(math::abs(D[i].imag()) <1.e-12 )
                            {
                                switch(settings_.spaceSolver){
                                case DIRECT:
                                    blockPrec->addOperator(i,i,gsComplexify<T>::make(gsSolverOp<sparseLLTfact>::make((D[i].real()*(*data_[level].Mx[slice]) + *data_[level].Kx[slice]))));
                                    blockEx->addOperator(i,i,blockPrec->getOperator(i,i));
                                    break;
                                case IETI:
                                {
                                    typename gsIETIAssembler<T>::Ptr ieti = memory::make_shared(new gsIETIAssembler<T>(*settings_.m_spatialAssemblers[spaceLevel]));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.reserve(data_[level].PatchMx[slice].size());
                                    for(size_t np=0; np<data_[level].PatchMx[slice].size();++np) {
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    D[i].real() * (*data_[level].PatchMx[slice][np])
                                                    + *data_[level].PatchKx[slice][np]);
                                        mats.push_back(give(mat));
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();
                                    typename IETIAdapter<T>::Ptr ptr = IETIAdapter<T>::make(ieti,settings_.precondIterations,settings_.precondTol,true,false,maxEigenvalue);
                                    blockPrec->addOperator(i,i,gsComplexify<T>::make(ptr));
                                    blockEx->addOperator(i,i,gsComplexify<T>::make(IETIAdapter<T>::make(ieti,100,1.e-10,false)));
                                    maxEigenvalue = ptr->getMaxEigenvalue();
                                }
                                    break;
                                default:
                                    GISMO_NO_IMPLEMENTATION;
                                };
                            }
                            else
                            {
                                for(int j=i+1;j<nt;++j)
                                {
                                    if(std::abs(D[i] - std::conj(D[j]))<1.e-12)
                                    {
                                        typename gsLinearOperator<T>::Ptr precond;

                                        switch(settings_.spaceSolver){
                                        case DIRECT:
                                            precond= gsSolverOp<sparseLLTfact>::make((gsSparseMatrix<T>(*data_[level].Kx[slice]+(D[i].real()+math::abs(D[i].imag()))**data_[level].Mx[slice])));
                                            break;
                                        case IETI:
                                        {
                                            typename gsIETIAssembler<T>::Ptr ieti = memory::make_shared(new gsIETIAssembler<T>(*settings_.m_spatialAssemblers[spaceLevel]));
                                            ieti->setOptions(settings_.IETIOptions);
                                            ieti->init();
                                            std::vector<gsSparseMatrix<T> > mats;
                                            mats.reserve(data_[level].PatchMx[slice].size());
                                            for(size_t np=0; np<data_[level].PatchMx[slice].size();++np) {
                                                gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                            (D[i].real() + math::abs(D[i].imag()))
                                                            * (*data_[level].PatchMx[slice][np])
                                                            + *data_[level].PatchKx[slice][np]);
                                                mats.push_back(give(mat));
                                            }
                                            ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                            ieti->assemble();
                                            precond = IETIAdapter<T>::make(ieti,settings_.precondIterations,settings_.precondTol,true,false,maxEigenvalue);
                                            maxEigenvalue =dynamic_cast<IETIAdapter<T>*>(precond.get())->getMaxEigenvalue();
                                        }
                                            break;
                                        };

                                        blockPrec->addOperator(i,i,ComplexSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],D[i],precond,settings_.approxIterates,settings_.tol,true,m_filename));
                                        blockEx->addOperator(i,i,ComplexSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],D[i],precond,100,1.e-10,false,m_filename));

                                        blockPrec->addOperator(j,j,ComplexSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],D[j],precond,settings_.approxIterates,settings_.tol,true,m_filename));
                                        blockEx->addOperator(j,j,ComplexSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],D[j],precond,100,1.e-10,false,m_filename));
                                        break;
                                    }


                                }
                            }
                        }
                        typename gsLinearOperator<std::complex<T> >::Ptr V = makePartialPivLUSolver(VinvTr);
                        precA[level][slice] = gsDeComplexify<T>::make(gsProductOp<std::complex<T> >::make(gsKroneckerOp<std::complex<T> >::make(V,gsIdentityOp<std::complex<T> >::make(nx)),blockPrec,gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(X),gsIdentityOp<std::complex<T> >::make(nx))));
                        LuOfA[level][slice]= gsDeComplexify<T>::make(gsProductOp<std::complex<T> >::make(gsKroneckerOp<std::complex<T> >::make(V,gsIdentityOp<std::complex<T> >::make(nx)),blockEx,gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(X),gsIdentityOp<std::complex<T> >::make(nx))));

                        /*
                    if(false)
                    {
                        gsSparseMatrix<T> A, temp;
                        temp=data[level].Kt[myTimeSlices[TimeLevels[level]][slice]]->kron(*data[level].Mx[slice],);
                        A=data[level].Mt[myTimeSlices[TimeLevels[level]][slice]]->kron(*data[level].Kx[slice]);
                        A+=temp;
                        LuOfA[level][slice] = gsSolverOp<sparseLUfact>::make(A);
                    }
                    */
                    }
                        break;
                    case COMPLEXSCHUR:
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

                            //(K + (alpha+|beta|)M)^{-1}
                            if(math::abs((*UT)(i,i).imag()) <1.e-12 )
                            {
                                switch(settings_.spaceSolver){
                                case DIRECT:
                                    APrec[i]=gsComplexify<T>::make(gsSolverOp<sparseLLTfact>::make(((*UT)(i,i).real()*(*data_[level].Mx[slice]) + *data_[level].Kx[slice])));
                                    AEx[i]=APrec[i];
                                    break;
                                case IETI:
                                {
                                    typename gsIETIAssembler<T>::Ptr ieti = memory::make_shared(new gsIETIAssembler<T>(*settings_.m_spatialAssemblers[spaceLevel]));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.reserve(data_[level].PatchMx[slice].size());
                                    for(size_t np=0; np<data_[level].PatchMx[slice].size();++np) {
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    (*UT)(i, i).real() * (*data_[level].PatchMx[slice][np])
                                                    + *data_[level].PatchKx[slice][np]);
                                        mats.push_back(give(mat));
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();

                                    typename IETIAdapter<T>::Ptr ptr = IETIAdapter<T>::make(ieti,settings_.precondIterations,settings_.precondTol,true,false,maxEigenvalue);
                                    APrec[i]=gsComplexify<T>::make(ptr);
                                    AEx[i]  = gsComplexify<T>::make(IETIAdapter<T>::make(ieti,100,1.e-10,false));
                                    maxEigenvalue  =ptr->getMaxEigenvalue();
                                }
                                    break;
                                default:
                                    GISMO_NO_IMPLEMENTATION;
                                };
                            }
                            else if(APrec[i]==NULL)
                            {
                                typename gsLinearOperator<T>::Ptr precond;
                                switch(settings_.spaceSolver){
                                case DIRECT:
                                    precond= gsSolverOp<sparseLLTfact>::make((gsSparseMatrix<T>(*data_[level].Kx[slice]+((*UT)(i,i).real()+math::abs((*UT)(i,i).imag()))**data_[level].Mx[slice])));
                                    break;
                                case IETI:
                                {
                                    typename gsIETIAssembler<T>::Ptr ieti = memory::make_shared(new gsIETIAssembler<T>(*settings_.m_spatialAssemblers[spaceLevel]));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.reserve(data_[level].PatchMx[slice].size());
                                    for(size_t np=0; np<data_[level].PatchMx[slice].size();++np) {
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    ((*UT)(i, i).real() + math::abs((*UT)(i, i).imag()))
                                                    * (*data_[level].PatchMx[slice][np]) + *data_[level].PatchKx[slice][np]);
                                        mats.push_back(give(mat));
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();

                                    precond = IETIAdapter<T>::make(ieti,settings_.precondIterations,settings_.precondTol,true, false,maxEigenvalue);
                                    maxEigenvalue  =dynamic_cast<IETIAdapter<T>*>(precond.get())->getMaxEigenvalue();
                                }
                                    break;
                                };

                                APrec[i]=ComplexSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],(*UT)(i,i),precond,settings_.approxIterates,settings_.tol,true,m_filename);
                                AEx[i]=ComplexSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],(*UT)(i,i),precond,100,1.e-10,false,m_filename);
                                for(int j=i+1;j<nt;++j)
                                {
                                    if(std::abs((*UT)(i,i) - std::conj((*UT)(j,j)))<1.e-12)
                                    {
                                        APrec[j]=ComplexSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],(*UT)(j,j),precond,settings_.approxIterates,settings_.tol,true,m_filename);
                                        AEx[j]=ComplexSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],(*UT)(j,j),precond,100,1.e-10,false,m_filename);
                                        break;
                                    }
                                }
                            }
                        }

                        precA[level][slice] = gsDeComplexify<T>::make( gsProductOp<std::complex<T> >::make(
                                                                                          gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(Pinv),gsIdentityOp<std::complex<T> >::make(nx)),
                                                                                          BlockTriangularSolver::make(APrec,data_[level].Mx[slice],UT),
                                                                                          gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(Q),gsIdentityOp<std::complex<T> >::make(nx))));

                        LuOfA[level][slice]=gsDeComplexify<T>::make(  gsProductOp<std::complex<T> >::make(
                                                                                        gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(Pinv),gsIdentityOp<std::complex<T> >::make(nx)),
                                                                                        BlockTriangularSolver::make(AEx,data_[level].Mx[slice],UT),
                                                                                        gsKroneckerOp<std::complex<T> >::make(makeMatrixOp(Q),gsIdentityOp<std::complex<T> >::make(nx))));

                        //     typename BlockTriangularSolver::Ptr solver =  BlockTriangularSolver::make(APrec,data_[level].Mx[slice],UT);
                        //       gsMatrix<std::complex<T> > result;
                        //      solver->toMatrix(result);

                        //        gsInfo<<"A: \n"<< A2.toDense()<<"\n\n";
                        //         gsInfo<<"A^{-1}: \n"<< A2.toDense().inverse()<<"\n\n";
                        //         gsInfo<<"pA: \n"<< result.inverse()<<"\n\n";
                        //          gsInfo<<"pA^{-1}: \n"<< result<<"\n\n";
                    }
                        break;
                    case REALSCHUR:
                    {
                        gsMatrix<T> Minv = Mt.inverse();
                        Eigen::RealSchur<Eigen::Matrix<T,Dynamic,Dynamic> > solv(Minv*Kt, true);
                        const gsMatrix<T >& TT = solv.matrixT();

                        gsMatrix<T> givens = gsMatrix<T>::Identity(nt,nt);

                        int nBlocks = 0;
                        for(int i=0;i<nt;++i)
                        {
                            if((i==nt-1 && math::abs(TT(i,i-1)) <1.e-12) || (i!=nt-1  && math::abs(TT(i+1,i)) <1.e-12) )
                                nBlocks++;
                            else
                            {
                                T t = (TT(i,i+1) + TT(i+1,i))/(TT(i,i)-TT(i+1,i+1));
                                //T t_s = math::min(math::abs(-t - math::sqrt(t*t + 1)),math::abs(-t + math::sqrt(t*t + 1)));
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

/*
                    gsInfo<<"Kt: \n"<<Kt<<"\n\n"<<"Mt: \n"<<Mt<<"\n\n";
                    gsInfo<<"Minv: \n"<<Minv<<"\n\n";
                    gsInfo<<"Q: \n"<<*Q<<"\n\n"<<"T: \n"<<*UT<<"\n\n"<<"Pinv: \n"<<*Pinv<<"\n\n";
                    gsInfo<<"Q*QT:\n "<<*Q*Q->transpose()<<"\n";
                    gsInfo<<"Givens:\n "<<givens<<"\n"<<std::flush;
                    gsInfo<<"Kt = PTQ^T:\n "<<Pinv->inverse()**UT*Q->transpose()<<"\n";

                    gsInfo<<"M*M^{-1} =M - P*Q^T:\n"<<Mt - Pinv->inverse()*Q->transpose()<<"\n\n";
                    gsInfo<<"PTQ^T  - K:\n"<<Pinv->inverse()**UT*Q->transpose() - Kt<<"\n\n";
//*/

                    /*
                    gsSparseMatrix<T > A,A2,D, temp;
                    temp=data_[level].Kt[myTimeSlices[TimeLevelsInp_[level]][slice]]->kron(*data_[level].Mx[slice]);
                    A=data_[level].Mt[myTimeSlices[TimeLevelsInp_[level]][slice]]->kron(*data_[level].Kx[slice]);
                    A+=temp;

                    temp=gsSparseMatrix<T>(UT->sparseView()).kron(*data_[level].Mx[slice]);
                    A2=gsSparseMatrix<T>(gsMatrix<T >::Identity(nt,nt).sparseView()).kron(*data_[level].Kx[slice]);
                    A2+=temp;
                    D=A2;

                    temp=gsSparseMatrix<T>(Q->transpose().sparseView()).kron(gsSparseMatrix<T>(gsMatrix<T>::Identity(nx,nx).sparseView()));
                    A2 = (A2*temp).eval();
                    temp=gsSparseMatrix<T>((Mt**Q).sparseView()).kron(gsSparseMatrix<T>(gsMatrix<T>::Identity(nx,nx).sparseView())); //Pinv->inverse().sparseView()
                    A2 = (temp*A2).eval();

                    gsInfo<<"A: \n"<<A.toDense()<<"\n\n";
                    gsInfo<<"A via decomp: \n"<<A2.toDense()<<"\n\n";
//*/
                        std::vector<typename gsLinearOperator<T >::Ptr> APrec(nBlocks);
                        std::vector<typename gsLinearOperator<T >::Ptr> AEx(nBlocks);
                        int idx = 0;
                        for(int i=0;i<nBlocks;++i)
                        {
                            //(K + (alpha+|beta|)M)^{-1}
                            if((i==nBlocks-1 && idx +1 ==nt) || (i!=nBlocks-1  && math::abs((*UT)(idx+1,idx)) <1.e-12) )
                            {
                                switch(settings_.spaceSolver){
                                case DIRECT:
                                    APrec[i]=gsSolverOp<sparseLLTfact>::make(((*UT)(idx,idx)*(*data_[level].Mx[slice]) + *data_[level].Kx[slice]));
                                    AEx[i]=APrec[i];
                                    break;
                                case IETI:
                                {
                                    typename gsIETIAssembler<T>::Ptr ieti = memory::make_shared(new gsIETIAssembler<T>(*settings_.m_spatialAssemblers[spaceLevel]));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.reserve(data_[level].PatchMx[slice].size());
                                    for(size_t np=0; np<data_[level].PatchMx[slice].size();++np) {
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    (*UT)(idx, idx) * (*data_[level].PatchMx[slice][np])
                                                    + *data_[level].PatchKx[slice][np]);
                                        mats.push_back(give(mat));
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();

                                    APrec[i] = IETIAdapter<T>::make(ieti,settings_.precondIterations,settings_.precondTol,true,false,maxEigenvalue);
                                    AEx[i] =IETIAdapter<T>::make(ieti,100,1.e-10,false);
                                    maxEigenvalue  =dynamic_cast<IETIAdapter<T>*>(APrec[i].get())->getMaxEigenvalue();
                                    /*
                                gsMatrix<T> mat, col;
                                mat.setZero(AEx[i]->cols(),AEx[i]->cols());
                                for(int c=0; c<AEx[i]->cols();++c)
                                {
                                    AEx[i]->apply(gsMatrix<T>::Identity(AEx[i]->cols(),AEx[i]->cols()).col(c),col);
                                    mat.col(c)=col;
                                }
                                for(size_t np=0; np<data_[level].PatchMx[slice].size();++np)
                                    gsInfo<<"Matrix on Patch: "<<np<<"\n"<<((*UT)(idx,idx)*(*data_[level].PatchMx[slice][np]) + *data_[level].PatchKx[slice][np]).toDense()<<"\n\n";

                                gsMatrix<T> orgA;
                                orgA = ((*UT)(idx,idx)*(*data_[level].Mx[slice]) + *data_[level].Kx[slice]).toDense();
                                gsInfo<<"original Matrix: \n"<<orgA<<"\n\n exact inverse: \n"<<orgA.inverse()<<"\n\n IETI approximation: \n"<<mat<<"\n\n IETIinv: \n"<<mat.inverse()<<"\n\n";
*/
                                }
                                    break;
                                default:
                                    GISMO_NO_IMPLEMENTATION;
                                };
                            }
                            else
                            {
                                //the abs is needed here, since the sign is not changed at this point
                                typename gsLinearOperator<T>::Ptr precond;
                                switch(settings_.spaceSolver){
                                case DIRECT:
                                    precond= gsSolverOp<sparseLLTfact>::make((gsSparseMatrix<T>(*data_[level].Kx[slice]+((*UT)(idx,idx) +sqrt(math::abs((*UT)(idx+1,idx)*(*UT)(idx,idx+1))))**data_[level].Mx[slice])));
                                    break;
                                case IETI:
                                {
                                    typename gsIETIAssembler<T>::Ptr ieti = memory::make_shared(new gsIETIAssembler<T>(*settings_.m_spatialAssemblers[spaceLevel]));
                                    ieti->setOptions(settings_.IETIOptions);
                                    ieti->init();
                                    std::vector<gsSparseMatrix<T> > mats;
                                    mats.reserve(data_[level].PatchMx[slice].size());
                                    for(size_t np=0; np<data_[level].PatchMx[slice].size();++np) {
                                        gsSparseMatrix<T> mat = gsSparseMatrix<T>(
                                                    ((*UT)(idx, idx) + sqrt(math::abs((*UT)(idx + 1, idx) * (*UT)(idx, idx + 1))))
                                                    * (*data_[level].PatchMx[slice][np]) + *data_[level].PatchKx[slice][np]);
                                        mats.push_back(give(mat));
                                    }
                                    ieti->giveAssembledMatrices(mats,gsMatrix<T>::Zero(ieti->getInfo().origSystemSize,1));
                                    ieti->assemble();

                                    precond = IETIAdapter<T>::make(ieti,settings_.precondIterations,settings_.precondTol,true,false,maxEigenvalue);
                                    maxEigenvalue  =dynamic_cast<IETIAdapter<T>*>(precond.get())->getMaxEigenvalue();
                                }
                                    break;
                                };

                                APrec[i]=RealSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],(*UT)(idx,idx),(*UT)(idx+1,idx),(*UT)(idx,idx+1),precond,settings_.approxIterates,settings_.tol,m_filename);
                                AEx[i]=RealSolver::make(data_[level].Kx[slice],data_[level].Mx[slice],(*UT)(idx,idx),(*UT)(idx+1,idx),(*UT)(idx,idx+1),precond,100,1.e-10,m_filename);
                                idx++;
                            }
                            idx++;
                        }

                        precA[level][slice] =gsProductOp<T>::make(   gsKroneckerOp<T>::make(makeMatrixOp(Pinv),gsIdentityOp<T>::make(nx)),
                                                                     QuasiBlockTriangularSolver::make(APrec,data_[level].Mx[slice],UT),
                                                                     gsKroneckerOp<T>::make(makeMatrixOp(Q),gsIdentityOp<T>::make(nx)));

                        LuOfA[level][slice]= gsProductOp<T>::make(  gsKroneckerOp<T>::make(makeMatrixOp(Pinv),gsIdentityOp<T>::make(nx)),
                                                                    QuasiBlockTriangularSolver::make(AEx,data_[level].Mx[slice],UT),
                                                                    gsKroneckerOp<T>::make(makeMatrixOp(Q),gsIdentityOp<T>::make(nx)));

                        /*
                    typename RealBlockTriangularSolver::Ptr solver =  RealBlockTriangularSolver::make(APrec,data_[level].Mx[slice],UT);


                    gsMatrix<T> result;
                    solver->toMatrix(result);

                    gsInfo<<"A: \n"<< D.toDense()<<"\n\n";
                    gsInfo<<"A^{-1}: \n"<< D.toDense().inverse()<<"\n\n";
                    gsInfo<<"pA: \n"<< result.inverse()<<"\n\n";
                    gsInfo<<"pA^{-1}: \n"<< result<<"\n\n";

                    gsInfo<<"error in solver: "<<(D.toDense().inverse()-result).norm()/(result.rows()*result.cols())<<std::endl;
*/
                    }
                        break;
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
        default:
            gsInfo<<"no other solver implemented, choose ("<<FULLSPACETIME<<") for direct solver or ("<<SPATIALWISE<<") for direct solver utilizing the tensor product structure.";
            GISMO_NO_IMPLEMENTATION;
            break;

        };

    }
}


//TODO: identicalSlice
template<typename T>
int gsHeatSTSlapOperator<T>::getNDofs(void) const { return ndofs_[0]; }

//TODO: identicalSlice
template<typename T>
bool gsHeatSTSlapOperator<T>::checkIfSpaceCoarseningIsAllowed(int timeLevel) const
{
    real_t criticalValue = 1.0;

    real_t tauL = (real_t)(math::exp2(NTimeLevels_-timeLevel-1))*tau_;

    if(tauL/(h_*h_) >= criticalValue) return true;
    else return false;
}

template<typename T>
void gsHeatSTSlapOperator<T>::mult(const Vector &u, Vector &f, int timeLevel, int timeStep, real_t sign) const
{
    //   f += sign*(A_[timeLevels_[timeLevel]].back()*u);
    gsMatrix<T> result;
    A[timeLevels_[timeLevel]][timeStep]->apply(u, result);
    // gsInfo<<"Vector before appl. A\n "<<u.transpose()<<"\n after appl. A\n "<<result.transpose()<<"\n";
    f +=sign*result;
}

//TODO: identicalSlice
template<typename T>
void gsHeatSTSlapOperator<T>::ApproximateSolve(Vector &u, const Vector &f, int timeLevel, int timeStep) const
{
    u.setZero();
    //gsInfo<<"before approximate solve: f.rows()= "<<f.rows()<<"\n";
    precA[timeLevels_[timeLevel]][timeStep]->apply(f,u);

}

template<typename T>
int gsHeatSTSlapOperator<T>::solveOnSTSlap(Vector &u, const Vector &f, int timeLevel, int timeStep) const
{
    int iterations = 0;
    LuOfA[timeLevels_[timeLevel]][timeStep]->apply(f,u);
    return iterations;
}

template<typename T>
void gsHeatSTSlapOperator<T>::calcInitialVector(const Vector &u, Vector &f, int timeLevel, int timeStep) const
{
    //GISMO_ASSERT(timeStep>0, "timeStep is 0, no B is defined here");
    Vector res;

    B[timeLevels_[timeLevel]][timeStep]->apply(u,res);
    // gsInfo<<"Vector before appl. B\n "<<u.transpose()<<"\n after appl. B\n "<<res.transpose()<<"\n";
    f+=res;
}

//Restrictions are time-slice wise
//Two timeslices are restricted to one
//On a timeslice the basis does not change (if all timeslices have the same basis).


//TODO: identicalSlice
template<typename T>
void gsHeatSTSlapOperator<T>::TimeRestriction1(const Vector & u_fine,
                                                     Vector & u_coarse,
                                                     int      /*coarseTimeLevel*/,
                                                     bool     setZero) const
{
    gsKroneckerOp<T> kron(makeMatrixOp(timeRestr1_),gsIdentityOp<T>::make((index_t)u_fine.rows()/timeRestr1_.cols()));
    if(setZero) u_coarse.setZero();
    // u_coarse += transferTime1[timeLevels_[coarseTimeLevel]]*u_fine;
    kron.apply(u_fine,m_temp);
    u_coarse += m_temp;
    /*    for(int i=0; i<m_; i++)
                                        {
                                            for(int k=0; k<blocksize_; k++)
                                            {
                                                int isk = k*m_+i;
                                                for(int l=0; l<blocksize_; l++)
                                                {
                                                    u_coarse[isk] += timeRestr1_(k,l)*u_fine[l*m_+i];
                                                }
                                            }
                                       }
                                   */
}
//TODO: identicalSlice
template<typename T>
void gsHeatSTSlapOperator<T>::TimeRestriction2(const Vector & u_fine,
                                                     Vector & u_coarse,
                                                     int      /*coarseTimeLevel*/,
                                                     bool     setZero) const
{
    gsKroneckerOp<T> kron(makeMatrixOp(timeRestr2_),gsIdentityOp<T>::make((index_t)u_fine.rows()/timeRestr2_.cols()));
    if(setZero) u_coarse.setZero();
    //u_coarse += transferTime2[timeLevels_[coarseTimeLevel]]*u_fine;
    kron.apply(u_fine,m_temp);
    u_coarse += m_temp;
    /*      for(int i=0; i<m_; i++)
                                        {
                                            for(int k=0; k<blocksize_; k++)
                                            {
                                                int isk = k*m_+i;
                                                for(int l=0; l<blocksize_; l++)
                                                {
                                                    u_coarse[isk] += timeRestr2_(k,l)*u_fine[l*m_+i];
                                                }
                                            }
                                        }
                                  */
}
//TODO: identicalSlice
template<typename T>
void gsHeatSTSlapOperator<T>::SpaceRestriction(const Vector & u_fine,
                                                     Vector & u_coarse,
                                                     int      coarseSpaceLevel,
                                                     bool     setZero) const
{
    if(setZero) u_coarse.setZero();
    //  gsInfo<<"ufine: "<<u_fine.rows()<<" - "<<u_fine.cols()<<";\tuCoarse: "<<u_coarse.rows()<<" - "<<u_coarse.cols()<<"\ttransferSpace: "<<transferSpace[coarseSpaceLevel]->rows()<<" - "<<transferSpace[coarseSpaceLevel]->cols()<<"\n"<<std::flush;
    transferSpace[coarseSpaceLevel]->apply(u_fine,m_temp);
    u_coarse += m_temp;
    //u_coarse += transferSpace[coarseSpaceLevel]*u_fine;
}

//TODO: identicalSlice
template<typename T>
void gsHeatSTSlapOperator<T>::TimeProlongation1(const Vector & u_coarse,
                                                      Vector & u_fine,
                                                      int      /*coarseTimeLevel*/,
                                                      bool     setZero) const
{
    gsKroneckerOp<T> kron(makeMatrixOp(timeRestr1_.transpose()),gsIdentityOp<T>::make((index_t)u_fine.rows()/timeRestr1_.cols()));
    if(setZero) u_fine.setZero();
    kron.apply(u_coarse,m_temp);
    u_fine += m_temp;

    //  u_fine+=  transferTime1[timeLevels_[coarseTimeLevel]].transpose()*u_coarse;
    /*       for(int i=0; i<m_; i++)
                                        {
                                            for(int k=0; k<blocksize_; k++)
                                            {
                                                int isk = k*m_+i;
                                                for(int l=0; l<blocksize_; l++)
                                                {
                                                    u_fine[isk] += timeRestr1_(l,k)*u_coarse[l*m_+i];
                                                }
                                            }
                                        }precA
                                 */
}

//TODO:identicalSlice
template<typename T>
void gsHeatSTSlapOperator<T>::TimeProlongation2(const Vector & u_coarse,
                                                      Vector & u_fine,
                                                      int      /*coarseTimeLevel*/,
                                                      bool     setZero) const
{
    gsKroneckerOp<T> kron(makeMatrixOp(timeRestr2_.transpose()),gsIdentityOp<T>::make((index_t)u_fine.rows()/timeRestr2_.cols()));
    if(setZero) u_fine.setZero();
    kron.apply(u_coarse,m_temp);
    u_fine += m_temp;

    // u_fine+=  transferTime2[timeLevels_[coarseTimeLevel]].transpose()*u_coarse;
    /*
                                        for(int i=0; i<m_; i++)
                                        {
                                            for(int k=0; k<blocksize_; k++)
                                            {
                                                int isk = k*m_+i;
                                                for(int l=0; l<blocksize_; l++)
                                                {
                                                    u_fine[isk] += timeRestr2_(l,k)*u_coarse[l*m_+i];
                                                }
                                            }
                                        }
                                        */
}

//TODO: identicalSlice
template<typename T>
void gsHeatSTSlapOperator<T>::SpaceProlongation(const Vector & u_coarse,
                                                      Vector & u_fine,
                                                      int      coarseSpaceLevel,
                                                      bool     setZero) const
{
    if(setZero) u_fine.setZero();
    transferTSpace[coarseSpaceLevel]->apply(u_coarse,m_temp);
    u_fine += m_temp;
    // u_fine += transferSpace[coarseSpaceLevel].transpose()*u_coarse;

}

template<typename T>
gsHeatSTSlapOperator<T>::ComplexSolver::ComplexSolver(const typename gsSparseMatrix<T>::Ptr   & K,
                                                      const typename gsSparseMatrix<T>::Ptr   & M,
                                                      const std::complex<T>                   & lambda,
                                                            typename gsLinearOperator<T>::Ptr   precond,
                                                            int                                 max_iter,
                                                            real_t                              tol,
                                                            bool                                /*prec*/,
                                                            std::string outputFilename) : m_K(K), m_M(M), m_lambda(lambda)
{
    m_filename = outputFilename;
    // m_MinRes = MinRes;
    m_max_iter = max_iter;
    m_tol = tol;

    m_inp.resize(2*cols(),1);
    m_out.resize(2*rows(),1);

    blockOP = gsBlockOp<T>::make(2,2);
    blockOP->addOperator(0,0, gsSumOp<T>::make(makeMatrixOp(K),gsScaledOp<T>::make(makeMatrixOp(M),lambda.real())));
    blockOP->addOperator(1,1, gsSumOp<T>::make(gsScaledOp<T>::make(makeMatrixOp(K),-1),gsScaledOp<T>::make(makeMatrixOp(M),-lambda.real())));
    blockOP->addOperator(0,1, gsScaledOp<T>::make(makeMatrixOp(M),lambda.imag()));
    blockOP->addOperator(1,0, gsScaledOp<T>::make(makeMatrixOp(M),lambda.imag()));

    blockPrec = gsBlockOp<T>::make(2,2);
    blockPrec->addOperator(0,0, precond);
    blockPrec->addOperator(1,1, precond);

    /*
    gsMatrix<T> big, pre;
    blockOP->toMatrix(big);
    blockPrec->toMatrix(pre);
    gsInfo<<"mat: \n"<<big<<std::endl;
    gsInfo<<"prec: \n"<<pre<<std::endl;
*/
    itSolver = memory::make_shared(new gsMinimalResidual<>(blockOP,blockPrec));
    itSolver->setMaxIterations(max_iter);
    itSolver->setTolerance(tol);
    m_print= true;
    //if(prec)
        itSolver->setInexactResidual(true);

}

template<typename T>
void gsHeatSTSlapOperator<T>::ComplexSolver::apply(const gsMatrix<std::complex<T> > & input, gsMatrix<std::complex<T> > & x) const
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
    /*
        else
    {
        blockPrec->apply(m_inp,m_out);
       // m_out = m_inp;
        m_out*=m_damping;
        for( int i=1; i<max_iter;++i )
        {
            m_out_old = m_out;
            blockOP->apply(m_out_old,m_temp);
            m_temp -= m_inp;
            blockPrec->apply(m_temp,m_out);
            //m_out = m_temp;
            m_out*=-m_damping;
            m_out +=  m_out_old;
        }

    }
*/
    x.real() = m_out.topRows(rows());
    x.imag() = -m_out.bottomRows(rows());
    gsInfo<<"\t Complex Solver numIT: "<<itSolver->iterations()<<"  - res: "<<itSolver->error()<<"\n"<<std::flush;
    if(m_print)
    {
        std::fstream out(m_filename.c_str(),std::ofstream::out | std::ofstream::app);
        out<<"Complex solver numIT: "<<itSolver->iterations()<<"  - res: "<<itSolver->error()<<"\n";
    }
}

template<typename T>
gsHeatSTSlapOperator<T>::RealSolver::RealSolver(const typename gsSparseMatrix<T>::Ptr& K, const  typename gsSparseMatrix<T>::Ptr& M, T a, T b, T c , typename gsLinearOperator<T>::Ptr precond, int max_iter, real_t tol, std::string outputFilename) : m_K(K), m_M(M), m_a(a), m_b(b), m_c(c)
{
    // m_MinRes = MinRes;
    m_max_iter = max_iter;
    m_tol = tol;
    m_filename = outputFilename;

    m_inp.resize(cols(),1);
    m_out.resize(rows(),1);

    blockOP = gsBlockOp<T>::make(2,2);
    /*
    blockOP->addOperator(0,0, gsSumOp<T>::make(makeMatrixOp(K),gsScaledOp<T>::make(makeMatrixOp(M),m_a)));
    blockOP->addOperator(1,1, gsSumOp<T>::make(gsScaledOp<T>::make(makeMatrixOp(K),-1),gsScaledOp<T>::make(makeMatrixOp(M),-m_a)));
    blockOP->addOperator(0,1, gsScaledOp<T>::make(makeMatrixOp(M),-m_c));
    blockOP->addOperator(1,0, gsScaledOp<T>::make(makeMatrixOp(M),m_b));

    blockPrec = gsBlockOp<T>::make(2,2);
    blockPrec->addOperator(0,0, precond);
    blockPrec->addOperator(1,1, precond);
*/

    blockOP->addOperator(0,0, gsSumOp<T>::make(gsScaledOp<T>::make(makeMatrixOp(K),math::abs(m_b)),gsScaledOp<T>::make(makeMatrixOp(M),math::abs(m_b)*m_a)));
    blockOP->addOperator(1,1, gsSumOp<T>::make(gsScaledOp<T>::make(makeMatrixOp(K),-math::abs(m_c)),gsScaledOp<T>::make(makeMatrixOp(M),-math::abs(m_c)*m_a)));
    blockOP->addOperator(0,1, gsScaledOp<T>::make(makeMatrixOp(M),-math::abs(m_b)*m_c));
    blockOP->addOperator(1,0, gsScaledOp<T>::make(makeMatrixOp(M),math::abs(m_c)*m_b));

    blockPrec = gsBlockOp<T>::make(2,2);
    blockPrec->addOperator(0,0,  gsScaledOp<T>::make(precond, T(1)/math::abs(m_b)));
    blockPrec->addOperator(1,1, gsScaledOp<T>::make(precond,T(1)/math::abs(m_c)));

    /*
    gsMatrix<T> big, pre;
    blockOP->toMatrix(big);
    blockPrec->toMatrix(pre);
    gsInfo<<"mat: \n"<<big<<std::endl;
    gsInfo<<"prec: \n"<<pre<<std::endl;
*/ //blockPrec
    //   itSolver = memory::make_shared(new gsGMRes<>(blockOP,blockPrec));
    itSolver = memory::make_shared(new gsMinimalResidual<real_t>(blockOP,blockPrec));
    
    itSolver->setMaxIterations(max_iter);
    itSolver->setTolerance(tol);
    itSolver->setInexactResidual(true);
    m_print= true;

}

template<typename T>
void gsHeatSTSlapOperator<T>::RealSolver::apply(const gsMatrix<T > & input, gsMatrix<T> & x) const
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
    gsInfo<<"\t Real Solver numIT: "<<itSolver->iterations()<<"  - res: "<<itSolver->error()<<"\n"<<std::flush;
    if(m_print)
    {
        std::fstream out(m_filename.c_str(),std::ofstream::out | std::ofstream::app);
        out<<"Real solver numIT: "<<itSolver->iterations()<<"  - res: "<<itSolver->error()<<"\n";
    }
}


template<typename T>
gsHeatSTSlapOperator<T>::BlockTriangularSolver::BlockTriangularSolver(std::vector<typename gsLinearOperator<std::complex<T> >::Ptr> A_inv, typename gsSparseMatrix<T>::Ptr B, typename gsMatrix<std::complex<T> >::Ptr UT) : m_Ainv(A_inv), m_B(makeMatrixOp(B)), m_T(UT)
{
    m_z.resize(A_inv.size(),gsMatrix<std::complex<T> >::Zero(B->rows(),1));
    m_rowVector.setZero(A_inv.size());
    for(size_t i=0; i<A_inv.size();++i )
        m_rowVector[i]=A_inv[i]->rows();
    m_colVector.setZero(1);

}
template<typename T>
gsHeatSTSlapOperator<T>::BlockTriangularSolver::BlockTriangularSolver(std::vector<typename gsLinearOperator<std::complex<T> >::Ptr> A_inv, typename gsLinearOperator<T>::Ptr B, typename gsMatrix<std::complex<T> >::Ptr UT) : m_Ainv(A_inv), m_B(B), m_T(UT)
{
    m_z.resize(A_inv.size(),gsMatrix<std::complex<T> >::Zero(B->rows(),1));
    m_rowVector.setZero(A_inv.size());
    for(size_t i=0; i<A_inv.size();++i )
        m_rowVector[i]=A_inv[i]->rows();
    m_colVector.setZero(1);

}

template<typename T>
void gsHeatSTSlapOperator<T>::BlockTriangularSolver::apply(const gsMatrix<std::complex<T> > & input, gsMatrix<std::complex<T> > & x) const
{
    m_colVector[0] = input.cols();
    x.setZero(rows(),input.cols());

    typename gsMatrix<std::complex<T> >::BlockView viewX = x.blockView(m_rowVector,m_colVector);
    for(int i= (int)m_Ainv.size()-1; i>=0;--i )
    {
        m_rhs = input.block(i*m_B->rows(),0,m_B->rows(),input.cols());
        if(i+1<(int)m_Ainv.size())
        {
            m_B->apply(viewX(i+1,0).real(),m_realTemp);
            m_z[i+1].real() = m_realTemp;
            m_B->apply(viewX(i+1,0).imag(),m_realTemp);
            m_z[i+1].imag() = m_realTemp;
        }

        for(size_t j=i+1;j<m_Ainv.size();++j)
            m_rhs-= (*m_T)(i,j)*m_z[j];

        m_Ainv[i]->apply(m_rhs,m_temp);
        viewX(i,0) = m_temp;

    }
}



template<typename T>
gsHeatSTSlapOperator<T>::QuasiBlockTriangularSolver::QuasiBlockTriangularSolver(std::vector<typename gsLinearOperator<T>::Ptr> A_inv, typename gsLinearOperator<T>::Ptr B, typename gsMatrix<T>::Ptr UT) : m_Ainv(A_inv), m_B(B), m_T(UT)
{
    m_z.resize(UT->cols(),gsMatrix<T>::Zero(B->rows(),1));
    m_rowVector.setZero(A_inv.size());
    for(size_t i=0; i<A_inv.size();++i )
        m_rowVector[i]=A_inv[i]->rows();
    m_colVector.setZero(1);

    GISMO_ASSERT(m_rowVector.sum() == rows(), "sizes do not match");

}

template<typename T>
gsHeatSTSlapOperator<T>::QuasiBlockTriangularSolver::QuasiBlockTriangularSolver(std::vector<typename gsLinearOperator<T>::Ptr> A_inv, typename gsSparseMatrix<T>::Ptr B, typename gsMatrix<T>::Ptr UT) : m_Ainv(A_inv), m_B(makeMatrixOp(B)), m_T(UT)
{
    m_z.resize(UT->cols(),gsMatrix<T>::Zero(B->rows(),1));
    m_rowVector.setZero(A_inv.size());
    for(size_t i=0; i<A_inv.size();++i )
        m_rowVector[i]=A_inv[i]->rows();
    m_colVector.setZero(1);

    GISMO_ASSERT(m_rowVector.sum() == rows(), "sizes do not match");
}

template<typename T>
void gsHeatSTSlapOperator<T>::QuasiBlockTriangularSolver::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    m_colVector[0] = input.cols();
    x.setZero(rows(),input.cols());

    typename gsMatrix<T>::BlockView viewX = x.blockView(m_rowVector,m_colVector);
    int start = input.rows();
    int entry = m_T->cols()-1;
    for(int i= (int)m_Ainv.size()-1; i>=0;--i )
    {
        start-=m_Ainv[i]->rows();
        m_rhs = input.block(start,0,m_Ainv[i]->rows(),input.cols());

        if((*m_B).rows()==m_Ainv[i]->rows())
        {
            for(index_t j=entry+1;j<(*m_T).cols();++j)
                m_rhs-= (*m_T)(entry,j)*m_z[j];

            m_Ainv[i]->apply(m_rhs,m_temp);
            viewX(i,0) = m_temp;
            entry--;

            if(entry >=0) m_B->apply(m_temp,m_z[entry+1]);

        }
        else
        {

            for(index_t j=entry+1;j<(*m_T).cols();++j)
            {
                m_rhs.topRows(m_B->rows())-= (*m_T)(entry-1,j)*m_z[j];
                m_rhs.bottomRows(m_B->rows())-= (*m_T)(entry,j)*m_z[j];
            }

            m_Ainv[i]->apply(m_rhs,m_temp);
            viewX(i,0) = m_temp;
            entry--;
            entry--;

            if(entry >=0)
            {
                m_B->apply(m_temp.topRows(m_B->rows()),m_z[entry+1]);
                m_B->apply(m_temp.bottomRows(m_B->rows()),m_z[entry+2]);
            }
        }
    }
}








} // namespace gismo

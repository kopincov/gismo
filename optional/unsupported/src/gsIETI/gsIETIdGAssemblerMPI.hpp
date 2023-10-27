/**  gsIETIdGAssemblerMPI.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on:  2015-05-07
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI
#include  <gsIETI/gsIETIdGAssemblerMPI.h>
#include  <gsAssembler/gsPoissonHeterogeneousAssembler.h>

namespace gismo {

template<class T>
void gsIETIdG_MPI_Base<T>::init_MPI(const gsIETIInfo& info,
const std::vector<std::pair<patchDof,patchDof> >& lagrangeTable,
std::vector<index_t>& freeElimLagr,
const std::vector<std::pair<patchDof,patchDof> >& elimExtraBasisConnection)
{
    gsIETI_MPI_Base<T>::init_MPI(info,lagrangeTable,freeElimLagr);

    m_ddof_Sbuffer.resize(m_comm.size());
    m_ddof_Rbuffer.resize(m_comm.size());
    m_ddof_request = new MPI_Request[infoMPI.nNeigbours];
    m_ddof_count = 0;
    std::vector<size_t> numbToS(info.numberPatches,0);
    std::vector<size_t> numbToR(info.numberPatches,0);
    m_ddof_connection.reserve(elimExtraBasisConnection.size());
    for( typename std::vector<std::pair<patchDof,patchDof> >::const_iterator  it = elimExtraBasisConnection.begin();it!=elimExtraBasisConnection.end();++it)
    {
        patchDof p1= (*it).first; //the original one
        patchDof p2 = (*it).second; //the extra Basis one

        if(Base::hasPatch(p1.first) || Base::hasPatch(p2.first))
            m_ddof_connection.push_back(*it);
    }

    for(typename std::vector<std::pair<patchDof,patchDof> >::const_iterator  it = m_ddof_connection.begin();it!=m_ddof_connection.end();++it)
    {
        patchDof p1= (*it).first; //the original one
        patchDof p2 = (*it).second; //the extra Basis one

        if(Base::hasPatch(p1.first) && !Base::hasPatch(p2.first))
            numbToS[p2.first]++;
        else if(!Base::hasPatch(p1.first) && Base::hasPatch(p2.first))
            numbToR[p1.first]++;
    }
    //  m_ddof_shiftS.resize(info.numberPatches);
    // m_ddof_shiftR.resize(info.numberPatches);
    std::vector<size_t> sumS(m_comm.size(),0),sumR(m_comm.size(),0);
    for(gsSortedVector<size_t>::const_iterator  np = infoMPI.patchNeigbour.begin();np!=infoMPI.patchNeigbour.end();++np)
    {
        // m_ddof_shiftS[np].reserve(math::pow(2,info.dim-1));
        // m_ddof_shiftR[neig].reserve(math::pow(2,info.dim-1));

        int p = m_patch2proc[*np];
        //  m_ddof_shiftS[neig].push_back(sumS);
        //   m_ddof_shiftR[neig].push_back(sumR);
        sumS[p]+=numbToS[*np];
        sumR[p]+=numbToR[*np];
    }

    for(unsigned p=0; p<infoMPI.nNeigbours;p++)
    {
        int proc = infoMPI.procNeigbour[p];
        m_ddof_Sbuffer[proc].setZero(sumS[proc],info.nRhs);
        m_ddof_Rbuffer[proc].setZero(sumR[proc],info.nRhs);

        //  gsDebug<<"proc: "<<m_comm.rank()<<" tries to sends to proc "<<proc<<"  with size " << sumS[proc]<<"\n";
        //  gsDebug<<"proc: "<<m_comm.rank()<<" awaits from proc "<<proc<<"  with size " << sumR[proc]<<"\n";
    }


}

template<class T>
void gsIETIdG_MPI_Base<T>::init_ddof()
{
    for(size_t p=0; p<infoMPI.nNeigbours;++p)
    {
        if(m_ddof_Rbuffer[infoMPI.procNeigbour[p]].size()!=0)
        {
            int proc = infoMPI.procNeigbour[p];
            MPI_Irecv(m_ddof_Rbuffer[proc].data(),m_ddof_Rbuffer[proc].size(),MPITraits<T>::getType(),proc,7,m_comm,&m_ddof_request[p]);
          //  gsDebug<<"proc: "<<m_comm.rank()<<"    wait from proc "<<proc<<"  with size: " << m_ddof_Rbuffer[proc].size()<<"\n";
            m_ddof_count++;
        }
        else
        {
            m_ddof_request[p] = MPI_REQUEST_NULL;
        }

    }
   // gsDebug<<"proc: "<<m_comm.rank()<<"  ddofcount: "<<m_ddof_count<<"\n";
}

template<class T>
void gsIETIdG_MPI_Base<T>::send_ddof()
{
    for(size_t p=0; p<infoMPI.nNeigbours;++p)
    {
        if(m_ddof_Sbuffer[infoMPI.procNeigbour[p]].size()!=0)
        {
            int proc = infoMPI.procNeigbour[p];
            MPI_Request req;
            MPI_Isend(m_ddof_Sbuffer[proc].data(),m_ddof_Sbuffer[proc].size(),MPITraits<T>::getType(),proc,7,m_comm,&req);
            // gsDebug<<"proc: "<<m_comm.rank()<<"    send to proc "<<proc<<"  with size: " << m_ddof_Sbuffer[proc].size()<<"\n";
            MPI_Request_free(&req);
        }
    }
}

template<class T>
bool gsIETIdG_MPI_Base<T>::finishOneDdof(int& p,gsMatrix<T>& ddof)
{
    int n;
    MPI_Waitany(infoMPI.nNeigbours,m_ddof_request,&n,MPI_STATUS_IGNORE);
    if(n==MPI_UNDEFINED || m_ddof_count==0)
        return false;

    p= infoMPI.procNeigbour[n];
    ddof= m_ddof_Rbuffer[p];
    m_ddof_count--;

   // gsDebug<<"proc: "<<m_comm.rank()<<"    got from proc "<<p<<"  with size: " << m_ddof_Rbuffer[p].size()<<"\n";

    return true;
    // np = infoMPI.patchNeigbour[n];

}

//-----------------------------------------------------------------------------------------------------

template<class T>
gsIETIdGAssemblerMPI<T>::gsIETIdGAssemblerMPI(gsAssembler<T>& assembler,gsSortedVector<size_t> myPatches, MPI_Comm comm)
    : gsIETIAssembler<T>(assembler),
      gsIETI_MPI_Base<T>(myPatches,comm),
      gsIETIdGAssembler<T>(assembler),
      gsIETIAssemblerMPI<T>(assembler,myPatches,comm),
      gsIETIdG_MPI_Base<T>(myPatches, comm)

{
    m_comm = gsMpiComm(comm);
}

template<class T>
gsOptionList gsIETIdGAssemblerMPI<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.update(dG_Base::defaultOptions(),gsOptionList::addIfUnknown);
    return opt;
}

template<class T>
void gsIETIdGAssemblerMPI<T>::setOptions(const gsOptionList &opt)
{
    gsIETIdGAssembler<T>::setOptions(opt);
    gsIETIAssemblerMPI<T>::setOptions(opt);
}


template<class T>
void gsIETIdGAssemblerMPI<T>::init()
{
    dG_Base::init();

    MPIdG_Base::init_MPI(info, m_lagrangeTable, m_freeElimLagr,m_elimExtraBasisConnection);

}




template<class T>
void gsIETIdGAssemblerMPI<T>::assemble(const gsMultiPatch<T> &curSol)
{

    Base::assembleInit();
    bool assembleRHS = !m_IETIoptions.opt.getSwitch("NoAssemblyOfRHS");
    MPI_Base::init_embT(info);
    MPI_Base::init_Spp(info);
    MPIdG_Base::init_ddof();

    std::vector<gsSparseMatrix<T> > matrices(info.numberPatches);
    std::vector<gsSparseMatrix<T> > matrices2(info.numberPatches);
    std::vector<gsMatrix<T> > rhs_loc(info.numberPatches);
   // std::vector<gsMatrix<T> > spp_loc(info.numberPatches);

    Base::reserveSpace(m_dirDofs,rhs_loc);

    Eigen::setNbThreads(1);

    gsStopwatch time, time2;
    size_t np;
    if(!m_IETIoptions.opt.getSwitch("ExtraTiming")) //The Fast One
    {
#pragma omp parallel
        {

            gsAssembler<T>* A = m_assembler->create();

#pragma omp for schedule(static, 1) nowait
            for (size_t i=0; i<m_patchIdx.size();++i)
            {
                np=m_patchIdx[i];
                Base::assembleLocal(A, matrices[np], rhs_loc[np],m_dirDofs[np],np,curSol);

                //gsInfo<<"K:\n\n"<<matrices[np].toDense()<<std::endl;

                //increase the size so that the extra dofs fit
                gsSparseMatrix<T>& matrix = matrices[np];
                matrix.uncompress();

                size_t newSize=0;
                for(size_t c=0; c<info.cDim;++c)
                    newSize+=m_locDofsMapper[np][c].freeSize();
                matrix.conservativeResize(newSize,newSize);

                //reserve for more dofs
                int nonZerosPerCol = m_assembler->numColNz();
                matrix.reserve(gsVector<index_t>::Constant(matrix.cols(), nonZerosPerCol));
                //   gsDebug<<"proc: "<<m_comm.rank()<<"  on patch "<<np<<":\n" <<m_dirDofs[np].transpose()<<"\n";
            }

            delete A;
#pragma omp barrier
            Base::printTime(time,"Time for assemling all patch local matrices: ");
            //fill up the remaining dirichlet dofs (only for the free-eliminated ones). Needed for assembleDgInterfaceContribution

#pragma omp single
            {

                prepareDDofs();
                MPIdG_Base::send_ddof();
                updateDDofs();
                Base::printTime(time,"Time for sending DDofs: ");

                //assemble DG contributions
                assembleDgInterfaceContribution(matrices,rhs_loc);

                // gsInfo<<"K:\n\n"<<matrices[0].toDense()<<std::endl;
                // gsInfo<<"K:\n\n"<<matrices[1].toDense()<<std::endl;
                Base::printTime(time,"Assemble dG Interface: ");

                for (size_t i=0; i<m_patchIdx.size();++i)
                    np=m_patchIdx[i], matrices[np].makeCompressed();
            }

#pragma omp for schedule(static, 1) nowait
            for (size_t i=0; i<m_patchIdx.size();++i)
            {
                np=m_patchIdx[i];
                // -- Reordering
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                {
                    matrices2[np] = matrices[np];
                    Base::makeReorderingPrimalRem(matrices2[np],np);
                }
                Base::makeReorderingBoundInt(matrices[np],np);


                // The main matrices for applying the system matrix
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                    Base::assembleKrrKrpKpp(matrices2[np],np);
                else
                {
                    dG_Base::assembleC(np);
                    Base::assembleLUofKC(matrices[np],np);
                }


                //assemble local Spp
                //Base::assembleSppLoc(MPI_Base::m_Spp_buffer[np], np);
                Base::assembleSppLoc(MPI_Base::m_Spp_sendBuffer, np);

            }

            MPI_Base::sendSPP();
            for (size_t i=0; i<m_patchIdx.size();++i)
            {
                np=m_patchIdx[i];
                //For Preconditioning and processing the rhs and solution ICO minimum energy
                Base::assembleKiiKibKbb(matrices[np], np);

                //assemble local rhs
                if(assembleRHS) Base::assembleRhs(rhs_loc, np);

            } //end-for
        }//end-parallel
    }//end if
    else // The slow one
    {

#pragma omp parallel
        {
            gsAssembler<T>* A = m_assembler->create();
#pragma omp for schedule(static, 1) nowait
            for (size_t i=0; i<m_patchIdx.size();++i)
            {
                np=m_patchIdx[i];
                Base::assembleLocal(A, matrices[np], rhs_loc[np],m_dirDofs[np], np,curSol);


                gsSparseMatrix<T>& matrix = matrices[np];
                matrix.uncompress();
                //increase the size so that the extra dofs fit

                size_t newSize=0;
                for(size_t c=0; c<info.cDim;++c)
                    newSize+=m_locDofsMapper[np][c].freeSize();
                matrix.conservativeResize(newSize,newSize);

                //reserve for more dofs
                int nonZerosPerCol = m_assembler->numColNz();
                matrix.reserve(gsVector<index_t>::Constant(matrix.cols(), nonZerosPerCol));

            }
            delete A;
        }

        Base::printTime(time,"Time for assemling all patch local matrices: ");

        //fill up the remaining dirichlet dofs (only for the free-eliminated ones)
        prepareDDofs();
        MPIdG_Base::send_ddof();
        updateDDofs();
        Base::printTime(time,"Time for sending DDofs: ");

        //assemble DG contributions
        assembleDgInterfaceContribution(matrices,rhs_loc);
        Base::printTime(time,"Assemble dG Interface: ");

        /*
        int max=0;
        for(unsigned np=0; np<info.numberPatches;np++)
        {
          //  gsInfo<<"the matrix starts here: \n"<<matrices[np].toDense()<<"\n\n";
            max=0;
            for(int i=0; i<matrices[np].cols();++i)
            {
                if(max<matrices[np].col(i).nonZeros())
                    max = matrices[np].col(i).nonZeros();

            }
            gsInfo<<"Max number of entries: "<<max <<"\n"<<"usage: "<<(real_t)max/m_assembler->options().numColNz(m_basis[0][0])<<"\n";
        }*/

#pragma omp parallel
        {
            // -- Reordering
#pragma omp for schedule(static, 1)
            for (size_t i=0; i<m_patchIdx.size();++i)
            {
                np=m_patchIdx[i];
                matrices[np].makeCompressed();
                if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
                {
                    matrices2[np] = matrices[np];
                    Base::makeReorderingPrimalRem(matrices2[np],np);
                }
                Base::makeReorderingBoundInt(matrices[np],np);
            }

            Base::printTime(time,"Time for reordering all patch local matrices: ");

            //------------------------------------------------------------------------------/
            //For Preconditioning and processing the rhs and solution ICO minimum energy



            // The main matrices for applying the system matrix
            if(m_IETIoptions.opt.getSwitch("NoMinimumEnergy"))
            {
#pragma omp for schedule(static, 1)
                for (size_t i=0; i<m_patchIdx.size();++i)
                    np=m_patchIdx[i],Base::assembleKrrKrpKpp(matrices2[np],np);

                Base::printTime(time,"Time for computing LU factorization of Krr : ");
            }
            else
            {
#pragma omp for schedule(static, 1)
                for (size_t i=0; i<m_patchIdx.size();++i)
                    np=m_patchIdx[i],dG_Base::assembleC(np);

                Base::printTime(time,"Time for computing C : ");

#pragma omp for schedule(static, 1)
                for (size_t i=0; i<m_patchIdx.size();++i)
                    np=m_patchIdx[i],Base::assembleLUofKC(matrices[np],np);

                Base::printTime(time,"Time for computing LU factorization of KC : ");
            }


#pragma omp for schedule(static, 1)
            for (size_t i=0; i<m_patchIdx.size();++i)
            {
               // np=m_patchIdx[i], Base::assembleSppLoc(MPI_Base::m_Spp_buffer[np], np);
                np=m_patchIdx[i], Base::assembleSppLoc(MPI_Base::m_Spp_sendBuffer, np);
            }
            MPI_Base::sendSPP();
            Base::printTime(time,"Time for assemling local Spp: ");

#pragma omp for schedule(static, 1)
            for (size_t i=0; i<m_patchIdx.size();++i)
                np=m_patchIdx[i],Base::assembleKiiKibKbb(matrices[np], np);

            Base::printTime(time,"Time for calculating LU of Kii : ");


#pragma omp for schedule(static, 1)
            for (size_t i=0; i<m_patchIdx.size();++i)
            {np=m_patchIdx[i]; if(assembleRHS) Base::assembleRhs(rhs_loc, np);}

            Base::printTime(time,"Time for assemling rhs: ");
        }
    }
    Base::printTime(time2,"Time for doing parallel assembling part: ");
    //------------Do Serial stuff--------------
    Eigen::setNbThreads(0);
    //1. right hand side
   if(assembleRHS) Base::assembleRhsFreeElim();

    time.restart();
    //2. Spp
    if(MPI_Base::isSppHolder())
    {
        MPI_Base::finalizeSpp(info, m_pDofsLoc2Glob,m_LU_Spp);
        printTime(time,"Time for assemling Spp: ");
      //  Base::finalizeEmbeddingTrans(m_rhs_p);
        // finalizeEmbeddingTrans(MPI_Base::m_embT_buffer,m_rhs_p);
      //  printTime(time,"Time for assemling final-rhs: ");
    }

    Base::printTime(time2,"Time for doing serial assembling part: ");
}

template<class T>
void gsIETIdGAssemblerMPI<T>::prepareDDofs()
{

    std::vector<size_t> ks(m_comm.size(),0);
    //fill up the remaining dirichlet dofs (only for the free-eliminated ones). Needed for assembleDgInterfaceContribution
    for(typename std::vector<std::pair<patchDof,patchDof> >::const_iterator it = MPIdG_Base::m_ddof_connection.begin(); it!=MPIdG_Base::m_ddof_connection.end() ;++it)
    {
        patchDof p1= (*it).first; //the original one
        patchDof p2 = (*it).second; //the extra Basis one

        if(MPI_Base::hasPatch(p1.first) && !MPI_Base::hasPatch(p2.first))
        {
            int c = getComp(p1.first,p1.second);
            int p= MPI_Base::m_patch2proc[p2.first];
            unsigned bIdx = m_locDofsMapper[p1.first][c].bindex(compCalcBack(p1.first,p1.second));
            MPIdG_Base::m_ddof_Sbuffer[p].row(ks[p]) = m_dirDofs[p1.first].row(bIdx);
            ks[p]++;
        }
        else if(MPI_Base::hasPatch(p1.first) && MPI_Base::hasPatch(p2.first))
        {
            int c = getComp(p2.first,p2.second); // p1 and p2 must have the same component

            unsigned bIdxExtra = m_locDofsMapper[p2.first][c].bindex(compCalcBack(p2.first,p2.second));
            unsigned bIdx = m_locDofsMapper[p1.first][c].bindex(compCalcBack(p1.first,p1.second));
            m_dirDofs[p2.first].row(bIdxExtra)= m_dirDofs[p1.first].row(bIdx);
        }


    }
    /*
    for(auto& p : infoMPI.procNeigbour)
        gsDebug<<"proc: "<<m_comm.rank()<<" sends to "<<p  <<"   "<<MPIdG_Base::m_ddof_Sbuffer[p].transpose()<<"\n";

    for(size_t np=0; np< info.numberPatches;np++)
        gsDebug<<"proc: "<<m_comm.rank()<<" has already on patch "<< np<<"   "<<m_dirDofs[np].transpose()<<"\n";
        */
}

template<class T>
void gsIETIdGAssemblerMPI<T>::updateDDofs()
{
    gsMatrix<T> ddof;

    int p;
    for(size_t neig=0; neig<infoMPI.nNeigbours;neig++)
    {
        if(!MPIdG_Base::finishOneDdof(p,ddof))
            continue;

        size_t k=0;
        for(typename std::vector<std::pair<patchDof,patchDof> >::const_iterator it = MPIdG_Base::m_ddof_connection.begin();it!=MPIdG_Base::m_ddof_connection.end();++it)
        {
            patchDof p1= (*it).first; //the original one
            patchDof p2 = (*it).second; //the extra Basis one
            if( MPI_Base::hasPatch(p2.first) && MPI_Base::m_patch2proc[p1.first]==p )
            {
                int c = getComp(p2.first,p2.second); // p1 and p2 must have the same component

                unsigned bIdxExtra = m_locDofsMapper[p2.first][c].bindex(compCalcBack(p2.first,p2.second));

              // gsDebug<<"proc: "<<m_comm.rank()<<"    want to insert index "<<k<<"  into index " << bIdxExtra<<"\n";
                m_dirDofs[p2.first].row(bIdxExtra)= ddof.row(k);
                k++;
            }

        }
    }

    /*
    for(size_t np=0; np< info.numberPatches;np++)
        gsDebug<<"proc: "<<m_comm.rank()<<" finishes on "<< np<<" with  "<<m_dirDofs[np].transpose()<<"\n";
*/
}

template<class T>
void gsIETIdGAssemblerMPI<T>::assembleDgInterfaceContribution(std::vector<gsSparseMatrix<T> >& matrices, std::vector<gsMatrix<T> > & rhs) const
{
    gsMatrix<index_t> *actives1,*actives2;
    gsMatrix<index_t> activesExtra1, activesExtra2;

    gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
    gsVector<T> quWeights;         // Mapped weights

    unsigned evFlags(0);     // Evaluation flags for the Geometry map
    gsQuadRule<T> QuRule;

    //TODO: make this stuff more generic.
    gsPoissonHeterogeneousAssembler<T>* pAss = static_cast< gsPoissonHeterogeneousAssembler<T>*>(m_assembler);
    GISMO_ASSERT(pAss !=NULL,
    "In order to use dG-IETI, the assembler needs to be derived from gsPoissonHeterogeneousAssembler<T>");

    int ii=0;
    for ( typename gsMultiPatch<T>::const_iiterator it = m_patches.iBegin(); it != m_patches.iEnd(); ++it )
    {
        if(dG_Base::m_iCoupling[ii++]!=iFace::dg)
            continue;

        const boundaryInterface & bi =
        ( m_basis[0][it->first() .patch].numElements(it->first() .side() ) <
        m_basis[0][it->second().patch].numElements(it->second().side() ) ?
        it->getInverse() : *it );

        const int patch1      = bi.first().patch;
        const int patch2      = bi.second().patch;

        if(!MPI_Base::hasPatch(patch1)&& !MPI_Base::hasPatch(patch2) )
            continue;

        const gsAffineFunction<T> interfaceMap(m_patches.getMapForInterface(bi));
        gsVisitorDg2<T,1> dg =  pAss->visitorDg(bi); //makes a copy, but class is very small.

        const gsBasis<T> & B1 = m_basis[0][patch1];// (!) unknown 0
        const gsBasis<T> & B2 = m_basis[0][patch2];

        const int bSize1      = B1.numElements( bi.first() .side() );
        const int bSize2      = B2.numElements( bi.second().side() );

        const int ratio = bSize1 / bSize2;
        GISMO_ASSERT(bSize1 >= bSize2 && bSize1%bSize2==0, "DG assumes nested interfaces. Got bSize1="<<bSize1<<", bSize2="<<bSize2<<"." );

        // Initialize
        dg.initialize(B1, B2, QuRule, evFlags);

        // Initialize geometry evaluators
        typename gsGeometryEvaluator<T>::uPtr geoEval1(getEvaluator(evFlags, m_patches[patch1]));
        typename gsGeometryEvaluator<T>::uPtr geoEval2(getEvaluator(evFlags, m_patches[patch2]));

        // Initialize domain element iterators
        typename gsBasis<T>::domainIter domIt1 = B1.makeDomainIterator( bi.first() .side() );
        typename gsBasis<T>::domainIter domIt2 = B2.makeDomainIterator( bi.second().side() );

        int count = 0;
        // iterate over all boundary grid cells on the "left"
        for (; domIt1->good(); domIt1->next() )
        {
            count++;
            // Compute the quadrature rule on both sides
            QuRule.mapTo( domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights);
            interfaceMap.eval_into(quNodes1,quNodes2);

            // Perform required evaluations on the quadrature nodes
            dg.evaluate(B1, *geoEval1, B2, *geoEval2, quNodes1, quNodes2);

            // Assemble on element
            if( MPI_Base::hasPatch(patch1) && MPI_Base::hasPatch(patch2))
                dg.assemble(*domIt1,*domIt2, *geoEval1, *geoEval2, quWeights);
            else
            {
                if(MPI_Base::hasPatch(patch2))
                {
                    dg.revert(true);
                    std::swap(geoEval1, geoEval2);
                    std::swap(domIt1, domIt2);
                }

                dg.assembleSingleSide(*domIt1,*domIt2, *geoEval1, *geoEval2, quWeights);

                if(MPI_Base::hasPatch(patch2))
                {
                    dg.revert(true);
                    std::swap(geoEval1, geoEval2);
                    std::swap(domIt1, domIt2);
                }
            }
            //extract the actives and prepare them for the IETI_locToGlob map
            dg.getActives(actives1, actives2);
            dG_Base::prepareActives(bi,*actives1,*actives2,activesExtra1,activesExtra2 );

            //do the map
            if(MPI_Base::hasPatch(patch1))
                dg.localToGlobalIETI(m_locDofsMapper[patch1][0],m_dirDofs[patch1], activesExtra1, matrices[patch1], rhs[patch1]);
            if(MPI_Base::hasPatch(patch2))
            {
                dg.revert();
                dg.localToGlobalIETI(m_locDofsMapper[patch2][0],m_dirDofs[patch2], activesExtra2, matrices[patch2], rhs[patch2]);
            }

            if ( count % ratio == 0 ) // next master element ?
            {
                domIt2->next();
            }

        }
    }

}

template<class T>
void gsIETIdGAssemblerMPI<T>::combineToCommonSolution(gsMatrix<T>& solVec) const
{
    m_comm.sum(solVec.data(), solVec.size()); //this is optional, usually there should be a class which stores a distrubuted function
}

} // namespace gismo

#endif

/**  gsIETISolverMPI.hpp
  
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofer
    Created on:  2014-12-02
*/

#include <gsCore/gsConfig.h>
#ifdef GISMO_WITH_MPI
#include <gsIETI/gsIETISolverMPI.h>


namespace gismo {


template<class T>
void gsIETISolverMPI<T>::applyMinEnergy(const gsMatrix<T> & inputP, const std::vector<gsMatrix<T> > & inputD, gsMatrix<T> &resultP, std::vector<gsMatrix<T> > &resultD,bool firstTime)const
{
    int nRhs = inputD[m_ass.m_patchIdx[0]].cols();
    Eigen::setNbThreads(1);
    gsMatrix<T> fP(info.dofTotalP,nRhs);
    fP.setZero();
   // m_ass.m_comm.barrier();

#pragma omp parallel
    {
        gsMatrix<T> rhs, sol;
#pragma omp for schedule(static, 1) nowait
        for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
        {
            size_t k=m_ass.m_patchIdx[i];

            rhs.setZero(info.dofsI[k]+info.dofsB[k]+info.dofsP[k], nRhs);
            sol.setZero(info.dofsI[k]+info.dofsB[k]+info.dofsP[k], nRhs);

            rhs.topRows(info.dofsB[k]) = inputD[k];
            m_ass.template solveKC<true,false>(k,rhs,sol);
            resultD[k]=sol.topRows(info.dofsB[k]);
            if(nonSymm && info.dofTotalP != 0)
                m_ass.applyPD(k,resultD[k],fP);
        }


        Eigen::setNbThreads(0);
    }

   // m_ass.m_comm.barrier();
    if(nonSymm && info.dofTotalP != 0)
        m_ass.m_comm.sum(fP.data(),fP.size()); //this could be a bit more efficient

    gsMatrix<T> iP_copy;
    //Finalize the sending of uPi and assemble uP
    if(!m_ass.hasFinalizedEmbT())
        m_ass.finalizeEmbeddingTrans(iP_copy,firstTime);
    else
        iP_copy = inputP;

    if(m_ass.isSppHolder())
    {
#pragma omp single nowait
        {
            iP_copy-=fP;
            if(info.dofTotalP != 0)
                m_ass.solveSpp(iP_copy,resultP);
        }

    }
 //   gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" after finalizing emb T with "<<resultP<<"\n";

  //  m_ass.initEmbedding(info, resultP);

}

template<class T>
void gsIETISolverMPI<T>::applyStdIETI(const gsMatrix<T> & inputP, const std::vector<gsMatrix<T> > & inputR, gsMatrix<T> &resultP, std::vector<gsMatrix<T> > &resultR,bool firstTime)const
{

    std::vector<gsMatrix<T> > vR(info.numberPatches);

    gsMatrix<T> vR_copy;
    gsMatrix<T> resP_temp, vP;

    int nRhs = inputR[m_ass.m_patchIdx[0]].cols();

    vP.setZero(info.dofTotalP,nRhs);

    if(m_ass.m_comm.rank()==0)
    {
        for(size_t np =0; np<info.numberPatches;++np)
        {
            if(m_ass.m_patch2proc[np]!=0)
            {
                m_vP_buff[np].setZero(info.dofsR[np],info.nRhs);
                MPI_Irecv(m_vP_buff[np].data(),info.dofsP[np]*info.nRhs,MPITraits<T>::getType(),m_ass.m_patch2proc[np],2,m_ass.m_comm,&m_requestVP[m_ass.m_patch2req[np]]);
            }
        }
    }


    //int indxR=0;
    for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
    {
        size_t k=m_ass.m_patchIdx[i];

        vR_copy = inputR[k];
        m_ass.solveKrr(k,vR_copy,vR[k]);
        m_vP_buff[k] = m_Krp[k]->transpose()*vR[k];
        if(m_ass.m_comm.rank()!=0)
        {
            MPI_Request req;
            MPI_Isend(m_vP_buff[k].data(),m_vP_buff[k].cols()*nRhs,MPITraits<T>::getType(),0,2,m_ass.m_comm,&req);
            MPI_Request_free(&req);
        }
    }



    //gsSparseMatrix<T> Kpr = m_Krp.transpose(); //TODO improve, e.g. K^Tx = (x^T K)^T
    gsMatrix<T> uP = gsMatrix<T>::Zero(info.dofTotalP,nRhs);
    m_ass.finalizeEmbeddingTrans(uP,firstTime);
    if(m_ass.isSppHolder())
    {
        MPI_Waitall(m_ass.m_req2patch.size(),m_requestVP,MPI_STATUSES_IGNORE);
        for(size_t np =0; np<info.numberPatches;++np)
            m_ass.assemblePrimal(m_vP_buff[np],np,vP);
        if(info.dofTotalP != 0)
        {
            gsMatrix<T> iP_copy = uP-vP;
            resultP.setZero(uP.rows(),uP.cols());
            m_ass.solveSpp(iP_copy,resultP);

        }
    }

    MPI_Bcast(resultP.data(),resultP.cols()*resultP.rows(),MPITraits<T>::getType(),0,m_ass.m_comm);

    gsMatrix<T> sol, temp;
    for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
    {
        size_t k=m_ass.m_patchIdx[i];

        resP_temp.setZero(info.dofsP[k],nRhs);
        m_ass.distributePrimal(resultP,k,resP_temp);

        if(info.dofsP[k]!=0)
        {
            temp.noalias() = *m_Krp[k]*resP_temp;
            m_ass.solveKrr(k,temp, sol);
            resultR[k]= vR[k] - sol;
        }
        else
            resultR[k]= vR[k];

    }

}

template<class T>
void gsIETISolverMPI<T>::calculateSolution(const gsMatrix<T> & lambda, gsMatrix<T> & resP, std::vector<gsMatrix<T> > & resR) const
{
    gsMatrix<T> fP,muP;
    std::vector<gsMatrix<T> >fR(info.numberPatches);
    std::vector<gsMatrix<T> > muR(info.numberPatches);

    m_ass.getRawRhs(fP,fR);

  //  if(gsMpi::worldRank()==0) gsInfo<<"primal rhs fP: "<<fP.transpose()<<"\n\n";


    muP.setZero(fP.rows(),fP.cols());
    m_B.applyTransTil(lambda, muP, muR);
   // gsInfo<<"fP: "<<fP.transpose()<<"\n";
 //   gsInfo<<"rank "<<m_ass.m_comm.rank()<<" finalizing embedding\n";
    //  if(m_ass.isSppHolder())
    m_ass.finalizeEmbeddingTrans(muP);
  //   if(gsMpi::worldRank()==0) gsInfo<<"primal muP: "<<muP.transpose()<<"\n\n";
   // gsInfo<<"muP: "<<muP.transpose()<<"\n";
  //  gsInfo<<"rank "<<m_ass.m_comm.rank()<<" finalizing embedding done\n";

    std::vector<gsMatrix<T> > diffR(info.numberPatches);

#pragma omp parallel for schedule(static, 1)
    for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
    {
        size_t k=m_ass.m_patchIdx[i];
        diffR[k]=fR[k]-muR[k];
    }
   // gsInfo<<"rank "<<m_ass.m_comm.rank()<<" before apply\n";
    if(m_ass.isSppHolder())
        fP-=muP;
    if(!m_ass.getOptions().opt.getSwitch("NoMinimumEnergy"))
        applyMinEnergy(fP,diffR,resP,resR);
    else
        applyStdIETI(fP,diffR,resP,resR);
   // gsInfo<<"rank "<<m_ass.m_comm.rank()<<" after apply\n";
}


//------------------------------------------------------------------------------------------------------------------------


template<class T>
gsIETIJumpOperatorMPI<T>::gsIETIJumpOperatorMPI(const gsIETIAssemblerMPI<T>& ass, IETIPrecondScaling::strategy sc)
: m_ass(ass), scaling(sc), m_info(ass.getInfo()), m_infoMPI(ass.getInfoMPI()){


    size_t sum= 0;
    for(size_t i=0; i< m_info.numberPatches;++i)
        sum+= m_info.dofsB[i];
    data.resize(sum);
    scalings.resize(m_info.numberPatches);
    size_t idx = 0;
    for(size_t i=0; i< m_info.numberPatches;++i)
    {
        scalings[i] = &(data[idx]);
        idx += m_info.dofsB[i];
    }

    if(scaling == IETIPrecondScaling::coefficient || scaling == IETIPrecondScaling::stiffnessModified)
    {
        m_coeffs.resize(m_info.numberPatches);
        if(scaling == IETIPrecondScaling::coefficient)
            setCoeffs();
        else
        {
            //use the average of the diagonal entries of the stiffness matrix.
            for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
            {
                size_t np=m_ass.m_patchIdx[i];
                T average=0;
                const gsSparseMatrix<T>& Kbb = *m_ass.getKbb(np);
                for(index_t k=0; k<m_info.dofsB[np];k++)
                    average+=Kbb(k,k);
                m_coeffs[np]=(average/m_info.dofsB[np]);
            }
        }
        /*
        for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
        {
            size_t np=m_ass.m_patchIdx[i];
            scalings[np].resize(m_info.dofsB[np]);
        }
*/
        for(typename std::map<index_t,std::set<patchDof> >::const_iterator it=m_ass.m_globalConnectionTable.begin(); it!=m_ass.m_globalConnectionTable.end(); ++it)
        {
            const std::set<patchDof>& set = it->second;

            T sumCoef = 0;
            for(typename std::set<patchDof>::const_iterator itS=set.begin();itS!=set.end();itS++)
            {
                int c = m_ass.getComp(itS->first,itS->second);
                if(m_ass.m_locDofsMapper[(*itS).first][c].is_free(m_ass.compCalcBack((*itS).first,(*itS).second)))
                    sumCoef+=m_coeffs[(*itS).first];
            }
            for(typename std::set<patchDof>::const_iterator itS=set.begin();itS!=set.end();itS++)
            {
                int c = m_ass.getComp(itS->first,itS->second);
                if(m_ass.m_locDofsMapper[(*itS).first][c].is_free(m_ass.compCalcBack((*itS).first,(*itS).second)))
                    scalings[itS->first][m_ass.m_glob2BoundInteriorIndex[itS->first][m_ass.m_locDofsMapper[itS->first][c].index(itS->second)]]= m_coeffs[itS->first]/sumCoef;
            }
        }

        /*
        for(size_t np=0; np<m_info.numberPatches; np++)
        {
            scalings[np].resize(m_info.dofsB[np]);
            for(unsigned k=0; k<m_info.dofsB[np];k++)
            {

                int c = m_ass.getComp(np,m_ass.m_boundDofs[np][k]);
                unsigned globalDof = m_ass.m_stdMapper[c].index(m_ass.compCalcBack(np,m_ass.m_boundDofs[np][k]),np);

                std::set<patchDof> set= m_ass.m_globalConnectionTable.find(globalDof)->second;
                T sumCoef = 0;
                for(typename std::set<patchDof>::const_iterator it=set.begin();it!=set.end();it++)
                {
                    if(m_ass.m_locDofsMapper[(*it).first][c].is_free(m_ass.compCalcBack((*it).first,(*it).second)))
                        sumCoef+=m_coeffs[(*it).first];
                }

                scalings[np][k]= m_coeffs[np]/sumCoef;

            }
        }
        */


    }
    else if(scaling == IETIPrecondScaling::multiplicity)
    {
        for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
        {
            size_t np=m_ass.m_patchIdx[i];

            // scalings[np].resize(m_info.dofsB[np]);
            for(index_t k=0; k<m_info.dofsB[np];++k)
            {
                int c = m_ass.getComp(np,m_ass.m_boundDofs[np][k]);
                unsigned globalDof = m_ass.m_stdMapper[c].index(m_ass.compCalcBack(np,m_ass.m_boundDofs[np][k]),np);
                scalings[np][k]= 1./(*m_ass.m_globalConnectionTable.find(globalDof)).second.size();
            }
        }
        m_ass.m_comm.sum(data.data(),data.size());
    }
    else if(scaling == IETIPrecondScaling::stiffness)
    {
        GISMO_ERROR("Stiffness scaling for MPI not implemented");
        const std::vector<typename gsSparseMatrix<T>::Ptr >& Kbb = m_ass.m_Kbb;
        for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
        {
            size_t np=m_ass.m_patchIdx[i];

            // scalings[np].resize(m_info.dofsB[np]);
            for(index_t k=0; k<m_info.dofsB[np];k++)
            {
                int c = m_ass.getComp(np,m_ass.m_boundDofs[np][k]);
                std::set<patchDof> set= m_ass.m_globalConnectionTable.find(m_ass.m_stdMapper[c].index(m_ass.compCalcBack(np,m_ass.m_boundDofs[np][k]),np))->second;
                T sumCoef = 0;
                for(typename std::set<patchDof>::const_iterator it=set.begin();it!=set.end();it++)
                {
                    unsigned patch = (*it).first;
                    int dof = m_ass.compCalcBack(patch,(*it).second);
                    if(m_ass.m_locDofsMapper[patch][c].is_free(dof))
                    {
                        index_t idx = m_ass.m_glob2BoundInteriorIndex[patch][m_ass.m_locDofsMapper[patch][c].index(dof)];
                        sumCoef+=(*Kbb[patch])(idx,idx);
                    }
                }
                scalings[np][k]= (real_t)(*Kbb[np])(k,k)/sumCoef;
            }
        }
        m_ass.m_comm.sum(data.data(),data.size());
    }



}

template<class T>
void gsIETIJumpOperatorMPI<T>::setCoeffs()
{
    m_coeffs.resize(m_info.numberPatches);

    for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
    {
        size_t np=m_ass.m_patchIdx[i];

        m_coeffs[np] = 1;

        const gsPdeWithCoeff<T>* ptr = dynamic_cast<const gsPdeWithCoeff<T>*>(&(m_ass.m_assembler->pde()));

        if( ptr )
            m_coeffs[np] = ptr->getCoeffForIETI(np);
        else
            gsWarn<<"setCoeffs: Assume homogeneous coefficient alpha==1.\n"
                "If the respective class has a member function .getCoeffForIETI(n), "
                "it must be derived from gsPdeWithCoeff<T>.\n\n";

        //if(m_ass.m_comm.rank()==0) gsInfo<<" coeff on patch "<<np<<": "<<m_coeffs[np]<<"\n";
    }
    m_ass.m_comm.sum(m_coeffs.data(),m_coeffs.size());
}

template<class T>
void gsIETIJumpOperatorMPI<T>::applyTransTil(const gsMatrix<T>& input,gsMatrix<T>&  uP, std::vector<gsMatrix<T> > &  u2) const
{
    std::vector<gsMatrix<T> > w(m_info.numberPatches);
    applyTrans(input, w);
   // gsInfo<<"after applyTrans "<<m_ass.m_comm.rank()<<"\n";

    m_ass.embeddingTrans(w,uP, u2);
    // gsInfo<<"after embeddingTrans "<<m_ass.m_comm.rank()<<"\n";
}

template<class T>
void gsIETIJumpOperatorMPI<T>::applyTil( gsMatrix<T>& xP, const std::vector<gsMatrix<T> > & x2, gsMatrix<T> & result) const
{
    std::vector<gsMatrix<T> > w(m_info.numberPatches);
    std::vector<gsMatrix<T> > xPis(m_ass.m_patchIdx.size());
    m_ass.initEmbedding(m_info,xP,xPis);
    m_ass.embedding(xPis,x2, w);
   // gsInfo<<"after embedding "<<m_ass.m_comm.rank()<<"\n";

    apply(w,result);
    //  gsInfo<<"after applyTil "<<m_ass.m_comm.rank()<<"\n";

}

template<class T>
void gsIETIJumpOperatorMPI<T>::apply(const std::vector<gsMatrix<T> >& w, gsMatrix<T> & result) const
{
    int rank = m_ass.m_comm.rank();
    int nRhs = w[m_ass.m_patchIdx[0]].cols();
    result.setZero(m_infoMPI.lagrangeMultReduce, nRhs);

    patchDof p1,p2;
    std::pair<patchDof,patchDof> entry;
    int sign;


    for(size_t i=0;i<m_ass.m_reducedLagrangeTable.size();i++){
        entry = m_ass.m_reducedLagrangeTable[i];
        p1 = entry.first;
        p2 = entry.second;

        p1.first > p2.first ? sign = 1 : sign = -1;

        int c = m_ass.getComp(p1.first,p1.second);//must be the same for both components

        int dof = m_ass.compCalcBack(p1.first,p1.second);
        if(m_ass.m_patch2proc[p1.first]==rank && m_ass.m_locDofsMapper[p1.first][c].is_free(dof))
        {
            T scal = getScaling(p2);
            int gIdx = m_ass.m_locDofsMapper[p1.first][c].index(dof);
            int a = m_ass.m_glob2BoundInteriorIndex[p1.first][gIdx];
            result.row(i) += sign*w[p1.first].row(a)*scal;
        }

        dof = m_ass.compCalcBack(p2.first,p2.second);
        if(m_ass.m_patch2proc[p2.first]==rank &&  m_ass.m_locDofsMapper[p2.first][c].is_free(dof))
        {
            T scal = getScaling(p1);
            int gIdx = m_ass.m_locDofsMapper[p2.first][c].index(dof);
            int b = m_ass.m_glob2BoundInteriorIndex[p2.first][gIdx];
            result.row(i) -= sign*w[p2.first].row(b)*scal;
        }

        //result.row(i)=w[p1.first].row(a) - w[p2.first].row(b);
    }
}

template<class T>
void gsIETIJumpOperatorMPI<T>::applyTrans(const gsMatrix<T>& input, std::vector<gsMatrix<T> >& w) const
{
    int rank = m_ass.m_comm.rank();
    int nRhs = input.cols();
    for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
    {
        size_t np=m_ass.m_patchIdx[i];
        w[np].setZero(m_info.dofsB[np],nRhs);
    }
    patchDof p1,p2;
    std::pair<patchDof,patchDof> entry;
    int sign;
    
    for(size_t i=0;i<m_ass.m_reducedLagrangeTable.size();i++)
    {
        entry = m_ass.m_reducedLagrangeTable[i];
        p1 = entry.first;
        p2 = entry.second;

        p1.first > p2.first ? sign = 1 : sign = -1;

        int c = m_ass.getComp(p1.first,p1.second); //must be the same for both components

        int dof = m_ass.compCalcBack(p1.first,p1.second);
        if(m_ass.m_patch2proc[p1.first]==rank &&  m_ass.m_locDofsMapper[p1.first][c].is_free(dof))
        {
            T scal = getScaling(p2);
            int a = m_ass.m_glob2BoundInteriorIndex[p1.first][m_ass.m_locDofsMapper[p1.first][c].index(dof)];
           // if(rank==0)gsInfo<<"rank: "<<rank<<" access "<<a<<" and patch"<<p1.first<<"with values: "<<sign<<" , "<<input(i,0)<<" , "<<scal<<"\n";
            w[p1.first].row(a) += sign*input.row(i)*scal ;

        }

        dof = m_ass.compCalcBack(p2.first,p2.second);
        if(m_ass.m_patch2proc[p2.first]==rank &&  m_ass.m_locDofsMapper[p2.first][c].is_free(dof))
        {
            T scal = getScaling(p1);
            int b=m_ass.m_glob2BoundInteriorIndex[p2.first][m_ass.m_locDofsMapper[p2.first][c].index(dof)];
           // if(rank==0)gsInfo<<"rank: "<<rank<<" access "<<b<<" and patch"<<p2.first<<"with values: "<<sign<<" , "<<input(i,0)<<" , "<<scal<<"\n";
            w[p2.first].row(b) -= sign*input.row(i)*scal;
        }


    }
}

template<class T>
inline T gsIETIJumpOperatorMPI<T>::getScaling(const patchDof& pd) const
{
    if(scaling != IETIPrecondScaling::none)
    {
        int dof = m_ass.compCalcBack(pd.first,pd.second);
        int c = m_ass.getComp(pd.first,pd.second);
        if(m_ass.m_locDofsMapper[pd.first][c].is_free(dof))
        {
            int idx = m_ass.m_glob2BoundInteriorIndex[pd.first][m_ass.m_locDofsMapper[pd.first][c].index(dof)];
            return scalings[pd.first][idx];
        }
        else
            return 1;
    }
    else
        return 1;
}

/*
template<class T>
void gsIETIJumpOperatorMPI<T>::calculateScaling(std::vector<gsMatrix<T> >& w) const
{
    if(scaling == IETIPrecondScaling::none)
        return;
    else
    {
        for(size_t np=0; np<m_info.numberPatches; np++)
            for(unsigned k=0; k<m_info.dofsB[np];k++)
                w[np].row(k)*=scalings[np][k];
    }
}
*/
template<class T>
void gsIETIJumpOperatorMPI<T>::calculateMatrixForm(gsMatrix<T>& B, gsMatrix<T>& BT)
{
    int nB=0;
    for(unsigned np = 0;np<m_info.numberPatches;np++)
        nB+=m_info.dofsB[np];

    int nL = m_info.lagrangeMult;

    B.setZero(nL,nB);

    std::vector< gsMatrix<T> >inp(m_info.numberPatches);
    int idx=0;
    for(size_t np = 0; np< m_info.numberPatches;np++)
    {
        inp[np].setZero(m_info.dofsB[np],nB);
        inp[np].block(0,idx,m_info.dofsB[np],m_info.dofsB[np]).setIdentity(m_info.dofsB[np],m_info.dofsB[np]);

        idx+=m_info.dofsB[np];
    }

    apply(inp,B);

    gsMatrix<T> id(nL,nL);
    id.setIdentity(nL,nL);
    std::vector<gsMatrix<T> > res(m_info.numberPatches);
    applyTrans(id,res);

    BT.setZero(nB,nL);
    idx=0;
    for(size_t np = 0; np< m_info.numberPatches;np++)
    {
        BT.block(idx,0,m_info.dofsB[np],nL)=res[np];
        idx+=m_info.dofsB[np];
    }
}

//------------------------------------------------------------------------------------------

} // namespace gismo

#endif

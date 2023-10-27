/**  gsIETISolver.hpp
  
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofer
    Created on:  2014-12-02
*/


#include <gsIETI/gsIETISolver.h>

#include <gsPde/gsPdeWithCoeff.h>

namespace gismo {


template<class T>
void gsIETISolver<T>::applyMinEnergy(const gsMatrix<T> & inputP, const std::vector<gsMatrix<T> > & inputD, gsMatrix<T> &resultP, std::vector<gsMatrix<T> > &resultD)const
{
    int nRhs = inputP.cols();
    Eigen::setNbThreads(1);
    gsMatrix<T> fP(inputP.rows(),inputP.cols());
    fP.setZero();
#pragma omp parallel
    {
        gsMatrix<T> rhs, sol;
        //    std::vector<gsMatrix<T> > resDP(info.numberPatches);

#pragma omp for schedule(static, 1) nowait //delete nowait
        for(size_t k = 0; k<info.numberPatches; k++)
        {
            rhs.setZero(info.dofsI[k]+info.dofsB[k]+info.dofsP[k], nRhs);
            sol.setZero(info.dofsI[k]+info.dofsB[k]+info.dofsP[k], nRhs);

            rhs.topRows(info.dofsB[k]) = inputD[k];
            m_ass. template solveKC<true,false>(k,rhs,sol);
            resultD[k]=sol.topRows(info.dofsB[k]);

            if(nonSymm)
                m_ass.applyPD(k,resultD[k],fP);
        }


        if(nonSymm)
        {
            #pragma omp barrier
        }


#pragma omp single nowait
        {
            if(info.dofTotalP != 0)
            {
                gsMatrix<T> iP_copy = inputP - fP;
                m_ass.solveSpp(iP_copy,resultP);
            }
        }



    }
    Eigen::setNbThreads(0);
}

template<class T>
void gsIETISolver<T>::applyStdIETI(const gsMatrix<T> & inputP, const std::vector<gsMatrix<T> > & inputR, gsMatrix<T> &resultP, std::vector<gsMatrix<T> > &resultR)const
{

    std::vector<gsMatrix<T> > vR(info.numberPatches);
    gsMatrix<T> vR_copy;
    gsMatrix<T> resP_temp, vP;

    int nRhs = inputP.cols();

    vP.setZero(info.dofTotalP,nRhs);

    //int indxR=0;
    for(size_t k = 0; k<info.numberPatches; k++)
    {
        vR_copy = inputR[k];
        m_ass.solveKrr(k,vR_copy,vR[k]);
        m_ass.assemblePrimal(m_Krp[k]->transpose()*vR[k],k,vP);
    }


    //gsSparseMatrix<T> Kpr = m_Krp.transpose(); //TODO improve, e.g. K^Tx = (x^T K)^T
    if(info.dofTotalP != 0)
    {
        gsMatrix<T> iP_copy = inputP-vP;
        m_ass.solveSpp(iP_copy,resultP);

    }

    gsMatrix<T> sol, temp;
    for(size_t k = 0; k<info.numberPatches; k++)
    {
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
void gsIETISolver<T>::calculateSolution(const gsMatrix<T> & lambda, gsMatrix<T> & resP, std::vector<gsMatrix<T> > & resR) const
{
    gsMatrix<T> fP;
    std::vector<gsMatrix<T> >fR(info.numberPatches);

    m_ass.getRawRhs(fP,fR);

    std::vector<gsMatrix<T> > muR(info.numberPatches);
    gsMatrix<T> muP(info.dofTotalP, info.nRhs);

    m_B.applyTransTil(lambda, muP, muR);

    std::vector<gsMatrix<T> > diffR(info.numberPatches);

#pragma omp parallel for schedule(static, 1)
    for(size_t k = 0; k<info.numberPatches; k++)
    {
        diffR[k]=fR[k]-muR[k];
    }

    if(!m_ass.getOptions().opt.getSwitch("NoMinimumEnergy"))
        applyMinEnergy(fP-muP,diffR,resP,resR);
    else
        applyStdIETI(fP-muP,diffR,resP,resR);

}
//--------------------------------------------------------------------------------

template<typename T>
gsInexactIETISolver<T>::gsInexactIETISolver(gsIETIAssembler<T>& assembler): m_B(assembler), m_ass(assembler)
{
    info=assembler.getInfo();

    m_Krp.resize(info.numberPatches);
    m_rhs.setZero(info.dofTotalP+info.dofTotalB+info.dofTotalI+info.lagrangeMult,info.nRhs);

    m_colBlocks.resize(1);
    m_rowBlocks.resize(2+2*info.numberPatches);
    for(size_t np=0; np<info.numberPatches;++np)
    {
        m_rowBlocks[2*np]=info.dofsB[np];
        m_rowBlocks[2*np+1]=info.dofsI[np];
    }
    m_rowBlocks[2*info.numberPatches]= info.dofTotalP;
    m_rowBlocks[2*info.numberPatches+1]=info.lagrangeMult;
    m_colBlocks[0] = info.nRhs;

    m_rowBlocksFull.resize(2+info.numberPatches);
    for(size_t np=0; np<info.numberPatches;++np)
        m_rowBlocksFull[np]=info.dofsB[np]+info.dofsI[np];
    m_rowBlocksFull[info.numberPatches]= info.dofTotalP;
    m_rowBlocksFull[info.numberPatches+1]=info.lagrangeMult;
}

template<typename T>
void gsInexactIETISolver<T>::init()
{
    typename gsMatrix<T>::BlockView blockVecRhs = m_rhs.blockView(m_rowBlocks,m_colBlocks);

    gsMatrix<T> rhs_p;
    std::vector<gsMatrix<T> > rhs_b(info.numberPatches);
    std::vector<gsMatrix<T> > rhs_i(info.numberPatches);
    m_ass.getRawRhs(rhs_p,rhs_b,rhs_i);

    if(!m_ass.nCoupledElimDofs())
        blockVecRhs(1+2*info.numberPatches)=- m_ass.getDirRhs(); //maybe -
    blockVecRhs(2*info.numberPatches) = rhs_p;

    for(size_t np=0; np<info.numberPatches;++np)
    {
        blockVecRhs(2*np) = rhs_b[np];
        blockVecRhs(2*np+1) = rhs_i[np];
    }

}

template<class T>
void gsInexactIETISolver<T>::applyMinEnergy(const gsMatrix<T>& input, gsMatrix<T>& output)const
{
    int nRhs = input.cols();
    Eigen::setNbThreads(1);
    gsVector<index_t,1> colBlocks;
    colBlocks<<input.cols();
    gsMatrix<T>& inp = const_cast<gsMatrix<T>&>(input);
    typename gsMatrix<T>::BlockView inp_block = inp.blockView(m_rowBlocksFull,colBlocks);
    typename gsMatrix<T>::BlockView out_block = output.blockView(m_rowBlocksFull,colBlocks);

#pragma omp parallel
    {
        gsMatrix<T> sol;
#pragma omp for schedule(static, 1) nowait
        for(size_t k = 0; k<info.numberPatches; k++)
        {
            sol.setZero(info.dofsB[k]+info.dofsI[k], nRhs);

            m_ass.applyK(k,inp_block(k),sol);

            //  m_ass.projectToDualSubspaceA(sol,k);
            //  resultB[k].conservativeResize(info.dofsB[k],resultB[k].cols());
            out_block(k) = sol;
        }

#pragma omp single nowait
        {
            gsMatrix<T> resultP;
            if(info.dofTotalP != 0)
                m_ass.applySpp(inp_block(info.numberPatches),resultP);
            out_block(info.numberPatches)= resultP;
        }
    }
    Eigen::setNbThreads(0);
}

template<class T>
void gsInexactIETISolver<T>::applyStdIETI(const gsMatrix<T>& input, gsMatrix<T>& output)const
{
    GISMO_UNUSED(input);
    GISMO_UNUSED(output);
    /*
    std::vector<gsMatrix<T> > vR(info.numberPatches);
    gsMatrix<T> vR_copy;
    gsMatrix<T> resP_temp, vP;

    int nRhs = inputP.cols();

    vP.setZero(info.dofTotalP,nRhs);

    //int indxR=0;
    for(size_t k = 0; k<info.numberPatches; k++)
    {
        vR_copy = inputB[k];
        m_ass.solveKrr(k,vR_copy,vR[k]);
        gsSparseMatrix<T> Kpr = m_Krp[k].transpose();
        m_ass.assemblePrimal(Kpr*vR[k],k,vP);
    }


    //gsSparseMatrix<T> Kpr = m_Krp.transpose(); //TODO improve, e.g. K^Tx = (x^T K)^T
    if(info.dofTotalP != 0)
    {
        gsMatrix<T> iP_copy = inputP-vP;
        m_ass.solveSpp(iP_copy,resultP);

    }

    gsMatrix<T> sol, temp;
    for(size_t k = 0; k<info.numberPatches; k++)
    {
        resP_temp.setZero(info.dofsP[k],nRhs);
        m_ass.distributePrimal(resultP,k,resP_temp);

        if(info.dofsP[k]!=0)
        {
            temp.noalias() = m_Krp[k]*resP_temp;
            m_ass.solveKrr(k,temp, sol);
            resultR[k]= vR[k] - sol;
        }
        else
            resultR[k]= vR[k];

    }
*/
}

template<class T>
void gsInexactIETISolver<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
{
    result.setZero(input.rows(),input.cols());
    gsVector<index_t> colBlocks(1);
    colBlocks<<input.cols();

    gsMatrix<T>& inp = const_cast<gsMatrix<T>&>(input);
    typename gsMatrix<T>::BlockView inp_block = inp.blockView(m_rowBlocks,colBlocks);
    typename gsMatrix<T>::BlockView out_block = result.blockView(m_rowBlocks,colBlocks);

    std::vector<gsMatrix<T> > u(info.numberPatches);
    for(size_t np=0; np<info.numberPatches;++np)
    {
        u[np].resize(info.dofsB[np]+info.dofsI[np],input.cols());
        u[np].topRows(info.dofsB[np])=inp_block(2*np);
        u[np].bottomRows(info.dofsI[np])=inp_block(2*np+1);
    }


    gsMatrix<T> x1P(info.dofTotalP,input.cols());
    std::vector<gsMatrix<T> > x1B(info.numberPatches);

    gsMatrix<T> lag_res(info.lagrangeMult,input.cols());

    m_B.applyTransTil(inp_block(2*info.numberPatches+1), x1P, x1B);
    m_B.applyTil(inp_block(2*info.numberPatches), u, lag_res);

    if(!m_ass.getOptions().opt.getSwitch("NoMinimumEnergy")){
        applyMinEnergy(input, result);
    }
    else{
        applyStdIETI(input, result);
    }

    out_block(2*info.numberPatches) += x1P;
    out_block(2*info.numberPatches+1)= lag_res;
    for(size_t np=0; np<info.numberPatches;++np)
        out_block(2*np) += x1B[np];

}

//------------------------------------------------------------------------------------------------------------------------


template<class T>
gsIETIJumpOperator<T>::gsIETIJumpOperator(const gsIETIAssembler<T>& ass, IETIPrecondScaling::strategy sc)
    : m_ass(ass), scaling(sc), m_info(ass.getInfo()){

    scalings.resize(m_info.numberPatches);

    if(scaling == IETIPrecondScaling::coefficient || scaling == IETIPrecondScaling::stiffnessModified)
    {
        m_coeffs.reserve(m_info.numberPatches);
        if(scaling == IETIPrecondScaling::coefficient)
            setCoeffs();
        else
        {
            //use the average of the diagonal entries of the stiffness matrix.
            for(unsigned np=0;np<m_info.numberPatches;np++)
            {
                T average=0;
                const gsSparseMatrix<T>& Kbb = *m_ass.getKbb(np);
                for(index_t k=0; k<m_info.dofsB[np];k++)
                    average+=Kbb(k,k);
                m_coeffs.push_back(average/m_info.dofsB[np]);
            }
        }
        for(size_t np=0; np<m_info.numberPatches; np++)
            scalings[np].resize(m_info.dofsB[np]);


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
        for(size_t np=0; np<m_info.numberPatches; np++)
        {
            scalings[np].resize(m_info.dofsB[np]);
            for(index_t k=0; k<m_info.dofsB[np];++k)
            {
                int c = m_ass.getComp(np,m_ass.m_boundDofs[np][k]);
                index_t globalDof = m_ass.m_stdMapper[c].index(m_ass.compCalcBack(np,m_ass.m_boundDofs[np][k]),np);
                scalings[np][k]= 1./(*m_ass.m_globalConnectionTable.find(globalDof)).second.size();
            }
        }
    }
    else if(scaling == IETIPrecondScaling::stiffness)
    {
        const std::vector<typename gsSparseMatrix<T>::Ptr >& Kbb = m_ass.m_Kbb;
        for(size_t np=0; np<m_info.numberPatches; np++)
        {
            scalings[np].resize(m_info.dofsB[np]);
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
    }
}

template<class T>
void gsIETIJumpOperator<T>::setCoeffs()
{
    unsigned nP = m_info.numberPatches;
    m_coeffs.reserve(m_info.numberPatches);

    for(unsigned np=0;np<nP;np++)
    {
        T coeff = 1;

        const gsPdeWithCoeff<T>* ptr = dynamic_cast<const gsPdeWithCoeff<T>*>(&(m_ass.m_assembler->pde()));

        if( ptr )
            coeff = ptr->getCoeffForIETI(np);
        else
            gsWarn<<"setCoeffs: Assume homogeneous coefficient alpha==1.\n"
                "If the respective class has a member function .getCoeffForIETI(n), "
                "it must be derived from gsPdeWithCoeff<T>.\n\n";

        if(coeff == 0)
            gsWarn<<"setCoeffs: The coefficient for the scaling is zero!\n";

        m_coeffs.push_back(coeff);
    }

    /*
    gsPoissonHeterogeneousAssembler<T>* assHet = dynamic_cast<gsPoissonHeterogeneousAssembler<T>*>(m_ass.m_assembler);
    if(assHet != NULL)
    {
        gsInfo<<"Using coefficients from the heterogenous Poisson with const coeff on each patch"<<std::endl;
        for(unsigned np=0;np<nP;np++)
            m_coeffs.push_back(assHet->getAlphaValue(np));
        //m_coeffs.push_back(m_ass.m_Kbb[np](0,0));
        return;
    }


    gsCDRAssembler<T>* ass0 = dynamic_cast<gsCDRAssembler<T>*>(m_ass.m_assembler);
    if(ass0 != NULL)
    {
        GISMO_NO_IMPLEMENTATION;
        return;
    }

    gsBiharmonicAssembler<T>* ass1 = dynamic_cast<gsBiharmonicAssembler<T>* >(m_ass.m_assembler);
    if(ass1 != NULL)
    {
        GISMO_NO_IMPLEMENTATION;
        return;
    }

    gsPoissonAssembler2<T>* ass2 = dynamic_cast<gsPoissonAssembler2<T>*>(m_ass.m_assembler);
    if(ass2 != NULL)
    {
        gsInfo<<"Using coefficients from the Poisson2 (alpha == 1)"<<std::endl;
        for(unsigned np=0;np<nP;np++)
            m_coeffs.push_back(1);
        return;
    }

    gsPoissonAssembler<T>* ass3 = dynamic_cast<gsPoissonAssembler<T>*>(m_ass.m_assembler);
    if(ass3 != NULL)
    {
        gsInfo<<"Using coefficients from the Poisson (alpha == 1)"<<std::endl;
        for(unsigned np=0;np<nP;np++)
            m_coeffs.push_back(1);
        return;
    }

    gsShellAssembler<T>* ass4 =dynamic_cast<gsShellAssembler<T>*>(m_ass.m_assembler);
    if(ass4 != NULL)
    {
        GISMO_NO_IMPLEMENTATION;
        return;
    }

    GISMO_ERROR("Your Assembler is not supported for this scaling method, choose a different one");
*/


}

template<class T>
void gsIETIJumpOperator<T>::applyTransTil(const gsMatrix<T>& input,gsMatrix<T>&  uP, std::vector<gsMatrix<T> > &  u2) const
{
    std::vector<gsMatrix<T> > w(m_info.numberPatches);
    applyTrans(input, w);

    m_ass.embeddingTrans(w,uP, u2);
}

template<class T>
void gsIETIJumpOperator<T>::applyTil(const gsMatrix<T>& xP, const std::vector<gsMatrix<T> > & x2, gsMatrix<T> & result) const
{
    std::vector<gsMatrix<T> > w(m_info.numberPatches);
    m_ass.embedding(xP,x2, w);

    apply(w,result);

}

template<class T>
void gsIETIJumpOperator<T>::apply(const std::vector<gsMatrix<T> >& w, gsMatrix<T> & result) const
{
    int nRhs = w.front().cols();
    result.setZero(m_info.lagrangeMult, nRhs);

    patchDof p1,p2;
    std::pair<patchDof,patchDof> entry;
    int sign;


    for(size_t i=0;i<m_ass.m_lagrangeTable.size();i++){
        entry = m_ass.m_lagrangeTable[i];
        p1 = entry.first;
        p2 = entry.second;

        p1.first > p2.first ? sign = 1 : sign = -1;

        int c = m_ass.getComp(p1.first,p1.second);//must be the same for both components

        int dof = m_ass.compCalcBack(p1.first,p1.second);
        if(m_ass.m_locDofsMapper[p1.first][c].is_free(dof))
        {
            T scal = getScaling(p2);
            int gIdx = m_ass.m_locDofsMapper[p1.first][c].index(dof);
            int a = m_ass.m_glob2BoundInteriorIndex[p1.first][gIdx];
            result.row(i) += sign*w[p1.first].row(a)*scal;
        }

        dof = m_ass.compCalcBack(p2.first,p2.second);
        if(m_ass.m_locDofsMapper[p2.first][c].is_free(dof))
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
void gsIETIJumpOperator<T>::applyTrans(const gsMatrix<T>& input, std::vector<gsMatrix<T> >& w) const
{

    int nRhs = input.cols();
    for(size_t np=0; np<m_info.numberPatches; np++)
    {
        w[np].setZero(m_info.dofsB[np],nRhs);
    }

    patchDof p1,p2;
    std::pair<patchDof,patchDof> entry;
    int sign;

    for(size_t i=0;i<m_ass.m_lagrangeTable.size();i++)
    {
        entry = m_ass.m_lagrangeTable[i];
        p1 = entry.first;
        p2 = entry.second;

        p1.first > p2.first ? sign = 1 : sign = -1;

        int c = m_ass.getComp(p1.first,p1.second); //must be the same for both components

        int dof = m_ass.compCalcBack(p1.first,p1.second);
        if(m_ass.m_locDofsMapper[p1.first][c].is_free(dof))
        {
            T scal = getScaling(p2);
            int a = m_ass.m_glob2BoundInteriorIndex[p1.first][m_ass.m_locDofsMapper[p1.first][c].index(dof)];
            w[p1.first].row(a) += sign*input.row(i)*scal ;
        }

        dof = m_ass.compCalcBack(p2.first,p2.second);
        if(m_ass.m_locDofsMapper[p2.first][c].is_free(dof))
        {
            T scal = getScaling(p1);
            int b=m_ass.m_glob2BoundInteriorIndex[p2.first][m_ass.m_locDofsMapper[p2.first][c].index(dof)];
            w[p2.first].row(b) -= sign*input.row(i)*scal;
        }


    }
}

template<class T>
inline T gsIETIJumpOperator<T>::getScaling(const patchDof& pd) const
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
void gsIETIJumpOperator<T>::calculateScaling(std::vector<gsMatrix<T> >& w) const
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
void gsIETIJumpOperator<T>::calculateMatrixForm(gsMatrix<T>& B, gsMatrix<T>& BT)
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

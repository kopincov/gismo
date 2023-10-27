/** @file gsIETIScaledDirichletMPI.h
  
    @brief
    
    
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofer
    Created on: 2015-01-19
*/

#pragma once
#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI

#include <gsIETI/gsDistributedOperator.h>
#include <gsIETI/gsIETIUtils.h>
#include <gsIETI/gsIETISolverMPI.h>
#include <gsIETI/gsIETIAssemblerMPI.h>


namespace gismo {

class gsIETIJumpOperatorMPI<class T>;

/**
 * @brief This class represents the preconditioner for the IETI method.
 * It supports several strategies for scaling:
 *      none , multiplicity, coefficient.
 * Note that the coefficient scaling is currently only available for poisson equation
 *
 *let H be the maximal diameter of a patch
 *  none: in the worst case leads to a condition number of the systemmatrix of O( (max(1+log(H_i/h_i))^3 )
 *  multiplicity, coefficient: leads to a condition number of the systemmatrix of O( (max(1+log(H_i/h_i))^2 )
 *
 */
template<class T>
class gsScaledDirichletPrecondMPI : public gsDistributedOperator<T>
{
public:
    /// Shared pointer for gsScaledDirichletPrecondMPI
    typedef memory::shared_ptr<gsScaledDirichletPrecondMPI> Ptr;

    /// Unique pointer for gsScaledDirichletPrecondMPI
    typedef memory::unique_ptr<gsScaledDirichletPrecondMPI> uPtr;
    
    /// Base class
    typedef memory::shared_ptr<gsDistributedOperator<T> > BasePtr;
    
    /// \brief shortcut for LU factorizations
    typedef typename gsIETIAssemblerMPI<T>::sparseLUfact sparseLUfact;

public:

    /**
     * @brief gsScaledDirichletPrecond The constructor for the preconditioner, requires the
     *  a gsIETIAssembler class and the strategy for scaling.
     * @param ass An instance of the gsIETIAssembler class
     * @param strat The desired method for scaling
     */
    gsScaledDirichletPrecondMPI(const gsIETIAssemblerMPI<T>& ass)
        : m_ass(ass),m_strat(ass.getOptions().scal) ,m_B(ass,ass.getOptions().scal), m_info(ass.getInfo()), m_infoMPI(ass.getInfoMPI())
    {
        //m_input_ac.setZero(m_infoMPI.lagrangeMultReduce,m_info.nRhs);
        m_Kib.resize(m_info.numberPatches);
        m_Kbi.resize(m_info.numberPatches);
        m_Kbb.resize(m_info.numberPatches);
        for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
        {
            size_t np=m_ass.m_patchIdx[i];
            m_Kib[np] = m_ass.getKib(np);
            m_Kbi[np] = m_ass.getKbi(np);
            m_Kbb[np] = m_ass.getKbb(np);
        }
    }
    
     static uPtr make( const gsIETIAssemblerMPI<T>& ass)
        { return memory::make_unique( new gsScaledDirichletPrecondMPI( ass ) ); }

    /**
     * @brief apply The main routine of the class, represents the application of the preconditioner
     * The algorithm is:
     * 1) apply the transposed jump operator (without embedding) input -> {w} on each patch
     * 2) scale the values of the boundary
     * 3) apply the Schur complement of the stiffness matrix (Kbb - Kbi*Kii^-1*Kib)
     * 4) scale the result again
     * 5) apply the jump operator  (without embedding) {w} on each patch -> {result}
     *
     * @param input lagrange multiplier (lagrangeMult x nRhs)
     * @param result the application of the preconditioner (lagrangeMult x nRhs)
     */
    void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
    {
        //gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" start apply\n";
        std::vector<gsMatrix<T> > w(m_info.numberPatches);

     //   m_input_ac= input;
        //Accumulate the distributed vector
       // gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" start accumulate\n"<<std::flush;
        accumulate(input,m_input_ac);

      //  gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" input_ac: "<<m_input_ac.transpose()<<"\n";

        //Apply transposed scaled jump operator
        //m_B.calculateScaling(w);

        m_B.applyTrans(m_input_ac,w);

      //  gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" w: "<<w[0].transpose()<<" , "<<w[1].transpose()<<"\n";

        //Apply S
        Eigen::setNbThreads(1);
#pragma omp parallel
        {
            gsMatrix<T> temp, wi;
#pragma omp for schedule(static, 1) nowait
            for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
            {
                size_t np=m_ass.m_patchIdx[i];
                temp.noalias() = *m_Kib[np]*w[np];
                m_ass.template solveKii<false>(np,temp,wi);

                w[np] = *m_Kbb[np]*w[np]-*m_Kbi[np]*wi;
            }
        }
        Eigen::setNbThreads(0);
        //Apply scaled jump operator
        //m_B.calculateScaling(w);

      //  gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" wNEW: "<<w[0].transpose()<<" , "<<w[1].transpose()<<"\n";
        //gsMatrix<T> temp;
        m_B.apply(w,result);

      //  gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" result_dist: "<<result.transpose()<<"\n";

        //Accumulate the distributed vector
      //

        accumulate(result,result);

     //   gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" result_acc: "<<result.transpose()<<"\n";

    }


    /// \brief return the number of rows
    index_t rows() const {return m_infoMPI.lagrangeMultReduce;}

    /// \brief return the number of columns
    index_t cols() const {return m_infoMPI.lagrangeMultReduce;}
/*
    const gsMatrix<T> & getAccumulated() const
    {
        return m_input_ac;
    }
*/

    void accumulate(const gsMatrix<T>& input, gsMatrix<T> & accumulated) const
    {
        m_ass.accumulate(input, accumulated);
     //   m_ass.accumulateCollective(input, accumulated);
    }


    inline void postAccumulate() const
    {
       m_ass.postAccumulate();
    }

    inline void startAccumulate(const  gsMatrix<T> & input) const
    {
       m_ass.startAccumulate(input);
    }

    inline void finishAccumulate(gsMatrix<T> & result) const
    {
        m_ass.finishAccumulate(result);
    }
    inline void buildGobalVector(const gsMatrix<T>& input, gsMatrix<T> & accumulated) const
    {
       m_ass.buildGlobalLagrangeMultiplier(input,accumulated);
    }
    virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
    {
        m_ass.distribute(input,distributed);
    }


    /**
     * @brief calculateMatrixForm The function provides the possibility to get the matrix representation of the apply function.
     * This is done by calling the apply function with all unit vector. The use is not recommended.
     * @param PF The matrix representation of the apply function. (lagrangeMult x lagrange Mult)
     */
    void calculateMatrixForm(gsMatrix<T>& PF)
    {
        PF.setZero(m_infoMPI.lagrangeMultReduce,m_infoMPI.lagrangeMultReduce);
        gsMatrix<T> I(m_infoMPI.lagrangeMultReduce,m_infoMPI.lagrangeMultReduce);
        I.setIdentity(m_infoMPI.lagrangeMultReduce,m_infoMPI.lagrangeMultReduce);

        apply(I,PF);
    }

    void calculateMatrixForm(std::vector<gsMatrix<T> >& S)
    {
        gsMatrix<> I,temp, wi, w_temp;
        for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
        {
            size_t np=m_ass.m_patchIdx[i];

            I.setIdentity(m_info.dofsB[np],m_info.dofsB[np] );

            temp.noalias() = *m_Kib[np]*I;
            m_ass.solveKii<false>(np,temp,wi);

            wi.transposeInPlace();

            w_temp.noalias() = wi*(*m_Kib[np]);
            w_temp.transposeInPlace();

            S.push_back(*m_Kbb[np]*I-w_temp);
        }

    }

    void calculateMatrixForm(gsMatrix<T>& B, gsMatrix<T>& BT)
    {
        int nB=0;
        for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
           nB+=m_info.dofsB[m_ass.m_patchIdx[i]];

        int nL = m_infoMPI.lagrangeMultReduce;

        B.setZero(nL,nB);

        std::vector< gsMatrix<T> >inp(m_info.numberPatches);
        int idx=0;
        for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
        {
            size_t np=m_ass.m_patchIdx[i];

            inp[np].setZero(m_info.dofsB[np],nB);
            inp[np].block(0,idx,m_info.dofsB[np],m_info.dofsB[np]).setIdentity(m_info.dofsB[np],m_info.dofsB[np]);

            idx+=m_info.dofsB[np];
        }

        m_B.apply(inp,B);

        gsMatrix<T> id(nL,nL);
        id.setIdentity(nL,nL);
        std::vector<gsMatrix<T> > res(m_info.numberPatches);
        m_B.applyTrans(id,res);

        BT.setZero(nB,nL);
        idx=0;
        for (size_t i=0; i<m_ass.m_patchIdx.size();++i)
        {
            size_t np=m_ass.m_patchIdx[i];

            BT.block(idx,0,m_info.dofsB[np],nL)=res[np];
            idx+=m_info.dofsB[np];
        }
    }
private:
    /// \brief reference of the gsIETIAssembler class
    const gsIETIAssemblerMPI<T>& m_ass;

    /// \brief strategy for scaling
    const IETIPrecondScaling::strategy& m_strat;

    /// \brief jump operator
    const gsIETIJumpOperatorMPI<T>  m_B;

    /// \brief info structure
    gsIETIInfo m_info;
    gsIETIInfoMPI m_infoMPI;

    /// \brief for each patch the matrix Kib
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kib;
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kbi;

    /// \brief for each patch the matrix Kbb
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kbb;

    mutable gsMatrix<T> m_input_ac;

   // mutable gsMatrix<T> m_input_ac;
};






} // namespace gismo
#endif

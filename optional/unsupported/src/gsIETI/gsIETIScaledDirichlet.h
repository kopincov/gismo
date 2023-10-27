/** @file gsIETIScaledDirichlet.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2015-01-19
*/


#pragma once

#include <gsSolver/gsLinearOperator.h>
#include <gsIETI/gsIETIUtils.h>
#include <gsIETI/gsIETISolver.h>
#include <gsIETI/gsIETIAssembler.h>


namespace gismo {


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
class gsScaledDirichletPrecond : public gsLinearOperator<>
{
public:
    /// Shared pointer for gsScaledDirichletPrecond
    typedef memory::shared_ptr<gsScaledDirichletPrecond> Ptr;

    /// Unique pointer for gsScaledDirichletPrecond
    typedef memory::unique_ptr<gsScaledDirichletPrecond> uPtr;

    /// Base class
    typedef memory::shared_ptr<gsLinearOperator<> > BasePtr;

    /// \brief shortcut for LU factorizations
    typedef typename gsIETIAssembler<T>::sparseLUfact sparseLUfact;

public:

    /**
     * @brief gsScaledDirichletPrecond The constructor for the preconditioner, requires the
     *  a gsIETIAssembler class and the strategy for scaling.
     * @param ass An instance of the gsIETIAssembler class
     * @param strat The desired method for scaling
     */
    gsScaledDirichletPrecond(const gsIETIAssembler<T>& ass)
        : m_ass(ass),m_strat(ass.getOptions().scal) ,m_B(ass,ass.getOptions().scal), m_info(ass.getInfo())
    {
        m_Kib.resize(m_info.numberPatches);
        m_Kbi.resize(m_info.numberPatches);
        m_Kbb.resize(m_info.numberPatches);
        for(size_t np = 0; np<m_info.numberPatches;np++)
        {
            m_Kib[np] = m_ass.getKib(np);
            m_Kbi[np] = m_ass.getKbi(np);
            m_Kbb[np] = m_ass.getKbb(np);
        }

    }

    static uPtr make( const gsIETIAssembler<T>& ass )
    { return memory::make_unique( new gsScaledDirichletPrecond( ass) ); }

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

        std::vector<gsMatrix<T> > w(m_info.numberPatches);


        //Apply transposed scaled jump operator
        //m_B.calculateScaling(w);
        m_B.applyTrans(input,w);


        //Apply S
        Eigen::setNbThreads(1);
#pragma omp parallel
        {
            gsMatrix<T> temp, wi, w_temp;
#pragma omp for schedule(static, 1) nowait
            for(int np =0; np<static_cast<int>(m_info.numberPatches);np++)
            {
                temp.noalias() = *m_Kib[np]*w[np];
                m_ass.template solveKii<false>(np,temp,wi);
               // wi.transposeInPlace();
               // w_temp.noalias() = wi*m_Kib[np];
              //  w_temp.transposeInPlace();
                w_temp = *m_Kbi[np]*wi;

                w[np] = *m_Kbb[np]*w[np]-w_temp;
            }
        }
        Eigen::setNbThreads(0);
        //Apply scaled jump operator
        //m_B.calculateScaling(w);
        m_B.apply(w,result);

    }


    /// \brief return the number of rows
    index_t rows() const {return m_info.lagrangeMult;}

    /// \brief return the number of columns
    index_t cols() const {return m_info.lagrangeMult;}

    /**
     * @brief calculateMatrixForm The function provides the possibility to get the matrix representation of the apply function.
     * This is done by calling the apply function with all unit vector. The use is not recommended.
     * @param PF The matrix representation of the apply function. (lagrangeMult x lagrange Mult)
     */
    void calculateMatrixForm(gsMatrix<T>& PF)
    {
        PF.setZero(m_info.lagrangeMult,m_info.lagrangeMult);
        gsMatrix<T> I(m_info.lagrangeMult,m_info.lagrangeMult);
        I.setIdentity(m_info.lagrangeMult,m_info.lagrangeMult);

        apply(I,PF);
    }

    void calculateMatrixForm(std::vector<gsMatrix<T> >& S)
    {
        gsMatrix<> I,temp, wi, w_temp;
        for(size_t np =0; np<m_info.numberPatches;np++)
        {
            I.setIdentity(m_info.dofsB[np],m_info.dofsB[np] );

            temp.noalias() = m_Kib[np]*I;
            GISMO_UNUSED(np);
            GISMO_UNUSED(temp);
            m_ass.solveKii<false>(np,temp,wi); //Clang 6.0.1 produces here a "warning: expression result unused:" for np and temp

          //  wi.transposeInPlace();
          //  w_temp.noalias() = wi*m_Kib[np];
          //  w_temp.transposeInPlace();

            w_temp = m_Kbi[np]*wi;

            S.push_back(m_Kbb[np]*I-w_temp);
        }

    }

    void calculateMatrixForm(gsMatrix<T>& B, gsMatrix<T>& BT)
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

        m_B.apply(inp,B);

        gsMatrix<T> id(nL,nL);
        id.setIdentity(nL,nL);
        std::vector<gsMatrix<T> > res(m_info.numberPatches);
        m_B.applyTrans(id,res);

        BT.setZero(nB,nL);
        idx=0;
        for(size_t np = 0; np< m_info.numberPatches;np++)
        {
            BT.block(idx,0,m_info.dofsB[np],nL)=res[np];
            idx+=m_info.dofsB[np];
        }
    }
private:
    /// \brief reference of the gsIETIAssembler class
    const gsIETIAssembler<T>& m_ass;

    /// \brief strategy for scaling
    const IETIPrecondScaling::strategy& m_strat;

    /// \brief jump operator
    const gsIETIJumpOperator<T>  m_B;

    /// \brief info structure
    gsIETIInfo m_info;

    /// \brief for each patch the matrix Kib
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kib;
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kbi;

    /// \brief for each patch the matrix Kbb
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Kbb;
};

namespace internal
{
template<typename T> class gsIETI_KinvWrapper;
}

template<class T>
class gsInexactIETIPrecond : public gsLinearOperator<>
{
public:
    /// Shared pointer for gsScaledDirichletPrecond
    typedef memory::shared_ptr<gsInexactIETIPrecond<T> > Ptr;

    /// Unique pointer for gsScaledDirichletPrecond
    typedef memory::unique_ptr<gsInexactIETIPrecond<T> > uPtr;

    /// Base class
    typedef memory::shared_ptr<gsLinearOperator<> > BasePtr;

    friend class internal::gsIETI_KinvWrapper<T>;

public:

    /**
     * @brief gsScaledDirichletPrecond The constructor for the preconditioner, requires the
     *  a gsIETIAssembler class and the strategy for scaling.
     * @param ass An instance of the gsIETIAssembler class
     * @param strat The desired method for scaling
     */
    gsInexactIETIPrecond(const gsIETIAssembler<T>& ass)
        : m_ass(ass), m_info(ass.getInfo()), m_Msd(gsScaledDirichletPrecond<T>::make(ass))
    {
        m_colBlocks.resize(1);
        m_rowBlocks.resize(2+m_info.numberPatches);
        for(size_t np=0; np<m_info.numberPatches;++np)
        {
            m_rowBlocks[np]=m_info.dofsB[np]+m_info.dofsI[np];
        }
        m_rowBlocks[m_info.numberPatches]= m_info.dofTotalP;
        m_rowBlocks[m_info.numberPatches+1]=m_info.lagrangeMult;

        m_colBlocks[0]=m_info.nRhs;

    }

    static Ptr make( const gsIETIAssembler<T>& ass )
    { return memory::make_shared( new gsInexactIETIPrecond( ass) ); }

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
        result.setZero(input.rows(),input.cols());
        gsVector<index_t> colBlocks(1);
        colBlocks<<input.cols();

        gsMatrix<T>& inp = const_cast<gsMatrix<T>&>(input);
        typename gsMatrix<T>::BlockView inp_block = inp.blockView(m_rowBlocks,colBlocks);
        typename gsMatrix<T>::BlockView out_block = result.blockView(m_rowBlocks,colBlocks);

        //apply preconditioner for Schur complement (Scaled dirichlet)
        gsMatrix<T> lag_out(m_info.lagrangeMult,input.cols());
        m_Msd->apply(inp_block(m_info.numberPatches+1),lag_out);
        out_block(m_info.numberPatches+1)=  lag_out;

        //---- Apply preconditioner for \tilde{K}

        //Block diagonal Preconditioner for K_\Delta\Delta
        gsMatrix<T> sol, rhs;
        for(size_t np=0; np<m_info.numberPatches;np++)
        {
            rhs.setZero(inp_block(np).rows()+m_info.dofsP[np],inp_block(np).cols());
            rhs.topRows(inp_block(np).rows()) = inp_block(np);
            m_ass. template solveKC<false,false>(np,rhs,sol);
            out_block(np) = sol.topRows(inp_block(np).rows());
        }

        gsMatrix<T> outP;
        m_ass.solveSpp(inp_block(m_info.numberPatches), outP);
        out_block(m_info.numberPatches)= outP;

    }


    /// \brief return the number of rows
    index_t rows() const {return  m_ass.systemSize();}

    /// \brief return the number of columns
    index_t cols() const {return  m_ass.systemSize();}

    /**
     * @brief calculateMatrixForm The function provides the possibility to get the matrix representation of the apply function.
     * This is done by calling the apply function with all unit vector. The use is not recommended.
     * @param PF The matrix representation of the apply function. (lagrangeMult x lagrange Mult)
     */
    void calculateMatrixForm(gsMatrix<T>& PF)
    {
        PF.setZero(rows(),cols());
        gsMatrix<T> I(rows(),cols());
        I.setIdentity(rows(),cols());

        apply(I,PF);
    }

    void calculateMatrixForm(std::vector<gsMatrix<T> >& S)
    {
        m_Msd->calculateMatrixForm(S);

    }

    void calculateMatrixForm(gsMatrix<T>& B, gsMatrix<T>& BT)
    {
        m_Msd->calculateMatrixForm(B,BT);
    }

    typename gsLinearOperator<T>::Ptr getPrecForK() const { return internal::gsIETI_KinvWrapper<T>::make(*this);}
    typename gsLinearOperator<T>::Ptr getPrecForS() const { return m_Msd;}
private:
    /// \brief reference of the gsIETIAssembler class
    const gsIETIAssembler<T>& m_ass;

    /// \brief info structure
    gsIETIInfo m_info;

    /// \brief Scaled dirichlet Preconditioner
    typename gsScaledDirichletPrecond<T>::Ptr m_Msd;


    gsVector<index_t> m_rowBlocks;
    gsVector<index_t> m_colBlocks;
};

namespace internal
{
template <typename T>
class gsIETI_KinvWrapper : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsScaledDirichletPrecond
    typedef memory::shared_ptr<gsIETI_KinvWrapper<T> > Ptr;

    /// Unique pointer for gsScaledDirichletPrecond
    typedef memory::unique_ptr<gsIETI_KinvWrapper<T> > uPtr;

    /// Base class
    typedef memory::shared_ptr<gsLinearOperator<> > BasePtr;

public:

    gsIETI_KinvWrapper(const gsInexactIETIPrecond<T>& iIETIP): m_iIETIP(iIETIP)
    {
        m_rowBlocks = m_iIETIP.m_rowBlocks.head(m_iIETIP.m_info.numberPatches+1);
    }

    static Ptr make(const gsInexactIETIPrecond<T>& iIETIP)
    { return memory::make_shared( new gsIETI_KinvWrapper<T>(iIETIP) );}

    void apply(const gsMatrix<T>& input, gsMatrix<T>& result) const
    {
        gsVector<index_t> colBlocks(1);
        colBlocks<<input.cols();

        gsMatrix<T>& inp = const_cast<gsMatrix<T>&>(input);
        typename gsMatrix<T>::BlockView inp_block = inp.blockView(m_rowBlocks,colBlocks);
        typename gsMatrix<T>::BlockView out_block = result.blockView(m_rowBlocks,colBlocks);

        //Block diagonal Preconditioner for K_\Delta\Delta
        gsMatrix<T> sol, rhs;
        for(size_t np=0; np<m_iIETIP.m_info.numberPatches;np++)
        {
            rhs.setZero(inp_block(np).rows()+m_iIETIP.m_info.dofsP[np],inp_block(np).cols());
            rhs.topRows(inp_block(np).rows()) = inp_block(np);

            m_iIETIP.m_ass. template solveKC<false,false>(np,rhs,sol);
            out_block(np) = sol.topRows(inp_block(np).rows());
        }

        gsMatrix<T> outP;
        m_iIETIP.m_ass.solveSpp(inp_block(m_iIETIP.m_info.numberPatches), outP);
        out_block(m_iIETIP.m_info.numberPatches)= outP;
    }

    /// \brief return the number of rows
    index_t rows() const {return  m_iIETIP.m_ass.systemSize()-m_iIETIP.m_ass.getInfo().lagrangeMult;}

    /// \brief return the number of columns
    index_t cols() const {return  m_iIETIP.m_ass.systemSize()-m_iIETIP.m_ass.getInfo().lagrangeMult;}

private:
    const gsInexactIETIPrecond<T>& m_iIETIP;
    gsVector<index_t> m_rowBlocks;
};
} // namespace internal


} // namespace gismo

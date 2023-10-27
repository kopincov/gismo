/** @file gsIETISolverMPI.h

    @brief


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2014-12-02
*/
#pragma once
#include <gsCore/gsConfig.h>

#ifdef GISMO_WITH_MPI

#include <gsIETI/gsIETIUtils.h>
#include <gsIETI/gsIETIAssemblerMPI.h>

#include <gsIETI/gsDistributedOperator.h>



namespace gismo {


/**
 *  This class represents the jump operator and calculates the jumps of a given function across the interfaces.
 *   It proived the routine for calling the the operator, its transposed and supports a scaled verion suited for precontioning.
 */
template<typename T>
class gsIETIJumpOperatorMPI
{
    typedef typename gsIETIAssemblerMPI<T>::patchDof patchDof; //patch + local number

    /// \brief  the corresponding IETI assembler
    const gsIETIAssemblerMPI<T>& m_ass;

    /// \brief  scaling strategy
    IETIPrecondScaling::strategy scaling;

    /// \brief reference to the info structure of gsIETIAssembler
    const gsIETIInfo & m_info;
    const gsIETIInfoMPI & m_infoMPI;

    /// \brief  contains for each interface dof on each patch its scaling
    std::vector<T*> scalings;
    std::vector<T> data;

    /// \brief  Stores for each patch the coefficient for scaling strategy
    std::vector<T> m_coeffs;

public:
    gsIETIJumpOperatorMPI(const gsIETIAssemblerMPI<T>& ass, IETIPrecondScaling::strategy sc = IETIPrecondScaling::none);

    /**
     * @brief applyTrans this function maps the lagrange multipliers (the jumps of the corresponding interface nodes)
     *   to the boundary values including a transformation to primal and remaining. This function combines applyBTrans(..) and gsIETIAssember.embeddingTrans
     * @param input the lagrange multipliers (lagrangeMult x nRhs)
     * @param uP the primal variables (totalPdofs x nRhs)
     * @param u2 the remaining variable for each patch (if isMinimumEnergy: dofsB[np] x nRhs ; else dofsR[np] x nRhs)
     */
    void applyTransTil(const gsMatrix<T>& input,gsMatrix<T>&  uP, std::vector<gsMatrix<T> > &  u2) const;

    /**
     * @brief apply this function maps the primal variables and the remaining to the jumnps across the interface, (= the lagrange multipliers).
     * This functions combines applyB and gsIETIAssember.embedding
     * @param xP the primal variables (totalPdofs x nRhs)
     * @param x2 the remaining variable for each patch (if isMinimumEnergy: dofsB[np] x nRhs ; else dofsR[np] x nRhs)
     * @param result the lagrange multipliers (lagrangeMult x nRhs)
     */
    void applyTil( gsMatrix<T>& xP, const std::vector<gsMatrix<T> > & x2, gsMatrix<T> & result) const;

    /**
     * @brief applyB this functions maps the values/coefficients on each patch to the jumps on the boundary
     * @param w for each patch it contains the inferface values/coefficients (dofsB x nRhs)
     * @param result the lagrange multipliers (lagrangeMult x nRhs)
     */
    void apply(const std::vector<gsMatrix<T> >& w, gsMatrix<T> & result) const;

    /**
     * @brief applyBTrans maps the lagrange multipliers (the jumps of the corresponding interface nodes)
     *   to the boundary values.
     * @param input the lagrange multipliers (lagrangeMult x nRhs)
     * @param w for each patch it contains the inferface values/coefficients (dofsB x nRhs)
     */
    void applyTrans(const gsMatrix<T>& input, std::vector<gsMatrix<T> >& w) const;


    //void calculateScaling(std::vector<gsMatrix<T> >& w) const;


    /**
     * @brief calculateMatrixForm calculates the matrix representation of applyB and applyBTrans included the scaling.
     * @param B the matrix representation for B
     * @param BT the matrix representation for B^T
     */
    void calculateMatrixForm(gsMatrix<T>& B, gsMatrix<T>& BT);

protected:

    /**
     * @brief calculateScaling given the opposite patch dof  (w.r.t. the lagrange multiplier), this function scales the values accoring to the set strategy for scaling. (see Dirichletpreconditioner).
     * If no preconditioning happens, all entries of scalings are 1, i.e. no scaling happends.
     * @param pd the opposite patch dof  (w.r.t. the lagrange multiplier)
     */
    inline T getScaling(const patchDof& pd) const;

    /**
     * @brief setCoeffs calculates the coefficients for scaling strategy. This factor is PDE dependent and may not be implemented.
     */
    void setCoeffs();
};





/**
 * \brief This class represents the application of the IETI system matrix and is
 * designed to fit as input for the gsConjugateGradient class. This means that this class is inherited from the abs It also returns the
 * corresponding right hand side.
 */
template<class T>
class gsIETISolverMPI : public gsDistributedOperator<T>
{
public:

    /// Shared pointer for gsIETISolverMPI
    typedef memory::shared_ptr<gsIETISolverMPI> Ptr;

    /// Unique pointer for gsIETISolverMPI
    typedef memory::unique_ptr<gsIETISolverMPI> uPtr;
    
    /// Base class
    typedef memory::shared_ptr<gsDistributedOperator<T> > BasePtr;
        
    /// \brief shortcut for LU-factorization
    typedef typename gsIETIAssemblerMPI<T>::sparseLUfact sparseLUfact;

    typedef  T Scalar;
    typedef real_t RealScalar;

protected:
    /**
     * @brief applyMinEnergy This funciton performs the application of the inverse Schur complement S
     * of the stiffness matrix via the minimal energy algorithm. It is part of the apply(..) function.
     *
     * @param inputP The primal part of the input argument. (dofsTotalP x nRhs)
     * @param inputD The dual part of the input argument. For each patch (dofsB x nRhs)
     * @param resultP The primal part of the result. (dofsTotalP x nRhs)
     * @param resultD The dual part of the result  For each patch (dofsB x nRhs)
     */
    void applyMinEnergy(const gsMatrix<T> & inputP, const std::vector<gsMatrix<T> > & inputD, gsMatrix<T> & resultP, std::vector<gsMatrix<T> > & resultD,bool firstTime=false) const;


    /**
     * @brief applyStdIETI his funciton performs the application of the inverse of the stiffness matrix K
     * via the the standard algorithm. It is part of the apply(..) function.
     *
     * @param inputP The primal part of the input argument. (dofsTotalP x nRhs)
     * @param inputR The remaining part of the input argument. For each patch (dofsR x nRhs)
     * @param resultP The primal part of the result. (dofsTotalP x nRhs)
     * @param resultR The remaining part of the result  For each patch (dofsB x nRhs)
     */
    void applyStdIETI(const gsMatrix<T> & inputP, const std::vector<gsMatrix<T> > & inputR, gsMatrix<T> & resultP, std::vector<gsMatrix<T> > & resultR,bool firstTime=false) const;


    /*
     *This functions are currently not used and needed.
     *
    void inline parseArray(const gsMatrix<T> & input, gsMatrix<T> &result1, std::vector<gsMatrix<T> >&result2) const
    {
        result1 = input.topRows(info.dofTotalP);
        parseStdVec(input.bottomRows(info.dofTotal-info.dofTotalP), result2);
    }

    void mergeArray(const gsMatrix<T> & x1,const std::vector<gsMatrix<T> > &x2, gsMatrix<T> &x) const
    {
        x.topRows(info.dofTotalP)= x1;
        mergeStdVec(x2, x.bottomRows(info.dofTotal-info.dofTotalP));
    }

    void inline parseStdVec(const gsMatrix<T>& input, std::vector<gsMatrix<T> >& result ) const
    {
        result.reserve(info.numberPatches);

        int indx = 0;
        for(index_t k=0; k<info.numberPatches;k++)
        {
            result.push_back(input.block(0,indx,info.dofs[k],info.nRhs));
            indx+= info.dofs[k];
        }
    }

    void inline mergeStdVec(const std::vector<gsMatrix<T> >& input, gsMatrix<T>& result )const
    {
        result.reserve(info.dofTotal);

        int indx = 0;
        for(index_t k=0; k<info.numberPatches;k++){
            result.block(0,indx,info.dofs[k],info.nRhs) = input[k];
            indx+= info.dofs[k];
        }
    }
    */


public:


    /**
     * @brief gsIETISolver Constructor of the IETISolver class.
     * @param assembler An instance of the gsIETIAssembler class. Based on this class gsIETISolver creates
     *  the jump operator B (gsIETIJumpOperator) and computes the appropriate rhs.
     */
    gsIETISolverMPI(gsIETIAssemblerMPI<T>& assembler): info(assembler.getInfo()), infoMPI(assembler.getInfoMPI()), m_B(assembler), m_ass(assembler)
    {
        m_rhs.resize(infoMPI.lagrangeMultReduce,info.nRhs);
        m_requestVP = new MPI_Request[m_ass.m_patch2req.size()];
    }

    static Ptr make(gsIETIAssemblerMPI<T>& assembler)
        { return memory::make_shared( new gsIETISolverMPI(assembler) ); }
    
    /**
     * @brief Destructor, does nothing. Contained pointers are freed in gsIETIAssembler
     */
    ~gsIETISolverMPI()
    {
         delete [] m_requestVP;
    }


    void init()
    {
        nonSymm = m_ass.getOptions().opt.getSwitch("NonSymmetric");
        m_rhs.setZero();
    //    if(m_ass.m_comm.rank() == 0 && nonSymm==true)
      //    gsInfo<<"Using non-symmetric version of IETI\n";

        gsMatrix<T> rhs_p, rhs_p_t;
        std::vector<gsMatrix<T> > rhs_r(info.numberPatches);
        std::vector<gsMatrix<T> > rhs_r_t(info.numberPatches);


        m_ass.getRawRhs(rhs_p,rhs_r);

        /*
        if(gsMpi::worldRank() == 0)
            gsInfo<<"fP: "<<rhs_p.transpose()<<"\n\n";
        for(size_t npi = 0; npi<info.numberPatches;++npi)
        {
            m_ass.m_comm.barrier();
            if(gsMpi::worldRank() == m_ass.m_patch2proc[npi])
                //gsInfo<<"patch: "<<np<<"  -- rhs_i: \t"<<m_rhs_i[np].transpose()<<"\n \n temp: \t"<<temp.transpose()<<"\n\n xI:\t"<<xI.transpose()<<"\n\n w:\t"<<w.transpose()<<"\n\n solVec:\t"<<solVec.transpose()<<"\n\n";
                gsInfo<<"fD "<<npi<<":\t"<<rhs_r[npi].transpose()<<"\n";
        }
        */

        // gsInfo<<"after getRawRhs "<<m_ass.m_comm.rank()<<"  rhs_p: "<<rhs_p.transpose()<<" rhs_r: "<<rhs_r[0].transpose()<<" , "<<rhs_r[1].transpose()<<"\n";

        if(!m_ass.getOptions().opt.getSwitch("NoMinimumEnergy")){
            //process rhs
            applyMinEnergy(rhs_p,rhs_r,rhs_p_t,rhs_r_t,true);
            // gsInfo<<"after apply() "<<m_ass.m_comm.rank()<<"  rhs_p_t: "<<rhs_p_t.transpose()<<" rhs_r_t: "<<rhs_r_t[0].transpose()<<" , "<<rhs_r_t[1].transpose()<<"\n";
        }
        else
        {
            m_vP_buff.resize(info.numberPatches);
            m_Krp.resize(info.numberPatches);

            for(size_t k = 0;k < info.numberPatches;k++)
                m_Krp[k] = m_ass.getKrp(k);

            applyStdIETI(rhs_p,rhs_r,rhs_p_t,rhs_r_t,true);
        }

        m_B.applyTil(rhs_p_t,rhs_r_t,m_rhs);
      //  gsInfo<<"after applyTil "<<m_ass.m_comm.rank()<<"\n";
        if(!m_ass.nCoupledElimDofs())
        {
        //   gsInfo<<"before handling extra ND "<<m_ass.m_comm.rank()<<"\n";
            m_rhs-= m_ass.getDirRhs();
          // gsInfo<<"after handling extra ND "<<m_ass.m_comm.rank()<<"\n";
        }
    }

    /**
     * @brief getRhs returns the rhs used for solving the system.
     * @return the corresponding rhs
     */
    gsMatrix<T>& getRhs()
    {
        return m_rhs;
    }

    gsMatrix<T> getRhs(index_t i)
    {
        GISMO_ASSERT(i<m_rhs.cols(), "You are trying to access a column which is not existent.");
        m_ass.setBufferCols(1); // resize the buffer to only one column, since we solve only for one rhs after the another
        return m_rhs.col(i);
    }

    /**
     * @brief apply The main routine of the class, performs the application of the system matrix. The application work in the following way:
     *
     * 1) B^T(input) -> u;
     * 2) applyMinEnergy(u) or applyStdIETI(u) -> v;
     * 3) B(v) -> result;
     *
     * The application of B includes:
     * 1) Embedding : inputP,inputD -> {input} on each patch
     * 2) Application of B: {input} on each patch -> result
     *
     * The application of B^T is vice verca
     * @param input The lagrange multipliers (lagrangeMult x nRhs)
     * @param result The application of F on the lagrange multipliers (lagrangeMult x Rhs)
     */
    void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const
    {
        gsMatrix<T> u1(info.dofTotalP,input.cols());
        std::vector<gsMatrix<T> > u2(info.numberPatches);

        gsMatrix<T> x1(info.dofTotalP,input.cols());
        std::vector<gsMatrix<T> > x2(info.numberPatches);

        m_B.applyTransTil(input, u1, u2);

      //  gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" before apply in solve\n ";
        if(!m_ass.getOptions().opt.getSwitch("NoMinimumEnergy")){
            applyMinEnergy(u1,u2, x1,x2);
        }
        else{
            applyStdIETI(u1, u2, x1, x2);
        }
       // gsInfo<<"rank: "<<m_ass.m_comm.rank()<<" after apply in solve\n ";
        m_B.applyTil(x1, x2, result);

    }

    /**
     * @brief rows returns the rows of the matrix. Inherited from base class
     * @return number of rows
     */
    index_t rows() const {return infoMPI.lagrangeMultReduce;}

    /**
     * @brief cols returns the columns of the matrix. Inherited from base class
     * @return number of columns
     */
    index_t cols() const {return infoMPI.lagrangeMultReduce;}


    /**
     * @brief calculateSolution After the solution of the system is found (e.g. via PCG), the lagrange multipliers are
     * used to calculate the solution in the form (primal, {dual} ). In order to receive a solution which can be plotted
     * use the function processSolution of gsIETIAssembler as with this data.
     * @param lambda The lagrange multipliers, i.e. solution where the PCG uses the apply function
     * @param resP The primal part of the solution
     * @param resR The dual/remaining part of the solution    /// \brief  Stores for each patch the coefficient for scaling strategy
    std::vector<T> m_coeffs;
     */
    void calculateSolution(const gsMatrix<T> & lambda, gsMatrix<T> & resP, std::vector<gsMatrix<T> > & resR) const;

    /**
     * @brief calculateSolution combines the functions calculateSolution(const gsMatrix<T> & lambda, gsMatrix<T> & resP, std::vector<gsMatrix<T> > & resR) and
     * gsIETIAssembler::processSolution(gsMatrix<T> &  resP,std::vector<gsMatrix<T> > & resR, gsMatrix<T> & solVec). Returns for given lagrange multipliers,
     * the solution of system in vector/matrix form. (like the output of the solver without using IETI)
     * Note that the member m_embT_buffer must be reset to its original size (nRhs)
     * @param lambda the lagrange multipliers
     * @param solVec the solution in matrix form
     */
    void calculateSolution(const gsMatrix<T> &lambda, gsMatrix<T>& solVec) const
    {
        m_ass.setBufferCols(lambda.cols());
        gsMatrix<T> resP;
        std::vector<gsMatrix<T> > resR(info.numberPatches);

        calculateSolution(lambda,resP,resR);
      //  if(gsMpi::worldRank()==0) gsInfo<<"primal Dofs: "<<resP.transpose()<<"\n\n";
        m_ass.processSolution(resP,resR,solVec);
    }

    /**
     * @brief calculateMatrixForm The function provides the possibility to get the matrix representation of the apply function.
     * This is done by calling the apply function with all unit vector. The use is not recommended.
     * @param F The matrix representation of the apply function. (lagrangeMult x lagrange Mult)
     */
    void calculateMatrixForm(gsMatrix<T>& F)
    {
        F.setZero(info.lagrangeMult,info.lagrangeMult);
        gsMatrix<T> I(info.lagrangeMult,info.lagrangeMult);
        I.setIdentity(info.lagrangeMult,info.lagrangeMult);

        apply(I,F);
    }

    void calculateMatrixForm(gsMatrix<T>& B, gsMatrix<T>& BT)
    {
        m_B.calculateMatrixForm(B,BT);
    }

    /*
    void calculateMatrixForm(std::vector<gsMatrix<T> >& S)
    {
        gsMatrix<T> IP,SP;
        std::vector<gsMatrix<T> > ID(info.numberPatches);
        std::vector<gsMatrix<T> > SD(info.numberPatches);


        std::vector<gsMatrix<T> >I(info.numberPatches);

        S.resize(info.numberPatches);
        for(unsigned np = 0; np< info.numberPatches;np++)
        {
            I[np].setIdentity(info.dofsB[np],info.dofsB[np]);
            S[np].setZero(info.dofsB[np],info.dofsB[np]);
         }

         m_ass.embeddingTrans(I,IP,ID);
         if(info.isMinimumEnergy){
             applyMinEnergy(IP,ID, SP,SD);
         }
         else{
             applyStdIETI(IP, ID, SP, SD);
         }

         m_ass.embedding(SP,SD,S);

    }
*/

    void accumulate(const gsMatrix<T>& input, gsMatrix<T> & accumulated) const
    {
        m_ass.accumulate(input, accumulated);
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


protected:
    /// \brief the info structure of gsIETIAssembler
    const gsIETIInfo & info;
    const gsIETIInfoMPI & infoMPI;

    /// \brief the jump operator
    gsIETIJumpOperatorMPI<T> m_B;

    /// \brief reference to the assembler class, if some informaton is required
    const gsIETIAssemblerMPI<T>& m_ass;

    /// \brief the right hand side of the system
    gsMatrix<T> m_rhs;

    /// \brief for each patch the matrix Krp
    std::vector<typename gsSparseMatrix<T>::Ptr > m_Krp;

    bool nonSymm;

protected:
     mutable std::vector<gsMatrix<T> > m_vP_buff;
     MPI_Request* m_requestVP;

};





} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIETISolverMPI.hpp)
#endif

#endif

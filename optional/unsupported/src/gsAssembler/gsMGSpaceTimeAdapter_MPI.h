/** @file gsMGSpaceTimeAdapter_MPI.h

    @brief Adapter Class for gsTimeParallelMultigrid class  considering the heat equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2017-11-23
*/


#pragma once

#include <gsAssembler/gsMGSpaceTimeAdapter.h>

#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpi.h>
#include <gsMpi/gsMpiComm.h>
#include <gsIETI/gsParallelOperator.h>
#include <gsIETI/gsParallelMinRes.h>
#include <gsIETI/gsDistributedBlockOp.h>


namespace gismo {

template<typename T>
class gsSpaceParallelHSTSlapOp : public gsHeatSTSlapOperator<T>
{
    typedef gsHeatSTSlapOperator<T> Base;
public:
    typedef typename  Base::localMatrix localMatrix;
    typedef typename Base::Vector Vector;

    /// Shared pointer for gsIdentityOp
    typedef memory::shared_ptr<gsSpaceParallelHSTSlapOp> Ptr;

    /// Unique pointer for gsIdentityOp
    typedef memory::unique_ptr<gsSpaceParallelHSTSlapOp> uPtr;

    typedef typename Base::STData STData;

public:
#if defined(GISMO_WITH_PARDISO)
    /*
    typedef typename gsSparseSolver<T>::PardisoLLT sparseLLTfact; //Threadsave
    typedef typename gsSparseSolver<T>::PardisoLU sparseLUfact; //Threadsave
    typedef typename gsSparseSolver<std::complex<T> >::PardisoLLT sparseLUfactC; //Threadsave
    */
    typedef typename gsSparseSolver<T>::SimplicialLDLT sparseLLTfact; //Threadsave
    typedef typename gsSparseSolver<T>::LU sparseLUfact; //Threadsave
    typedef typename gsSparseSolver<std::complex<T> >::SimplicialLDLT sparseLUfactC; //Threadsave
#elif defined(_OPENMP)
    typedef typename gsSparseSolver<T>::SimplicialLDLT sparseLLTfact; //Threadsave
    typedef typename gsSparseSolver<T>::LU sparseLUfact; //Threadsave
    typedef typename gsSparseSolver<std::complex<T> >::SimplicialLDLT sparseLUfactC; //Threadsave
    //#elif(defined(GISMO_WITH_SUPERLU))
    //    typedef typename gsSparseSolver<T>::SuperLU sparseLLTfact; //Not Threadsave
    //    typedef typename gsSparseSolver<T>::SuperLU sparseLUfact; //Not Threadsave
    //    typedef typename gsSparseSolver<std::complex<T> >::SuperLU sparseLUfactC; //Not Threadsave
#else
    typedef typename gsSparseSolver<T>::SimplicialLDLT sparseLLTfact; //Not Threadsave
    typedef typename gsSparseSolver<T>::LU sparseLUfact; //Threadsave
    typedef typename gsSparseSolver<std::complex<T> >::SimplicialLDLT sparseLUfactC; //Threadsave
#endif


/*
    class IETIApplicationMPI : public gsDistributedOperator<T>
    {
    private:
        class accumulator
        {
        public:
            accumulator(const typename gsScaledDirichletPrecondMPI<T>::Ptr& IETIPrec) : m_IETIPrec(IETIPrec) {}

            void accumulate(const gsMatrix<T> &input, gsMatrix<T> &x) const
            {
                m_IETIPrec->accumulate(input,x);
            }

        private:

            mutable typename gsScaledDirichletPrecondMPI<T>::Ptr m_IETIPrec;
        };

    public:
        /// Shared pointer for gsLinearOperator
        typedef memory::shared_ptr<IETIApplicationMPI > Ptr;

        /// Unique pointer for gsLinearOperator
        typedef memory::unique_ptr<IETIApplicationMPI > uPtr;


        IETIApplicationMPI(typename gsIETIAssemblerMPI<T>::Ptr IETIAss): m_IETIAss(IETIAss) , m_IETISolv(gsIETISolverMPI<T>::make(*m_IETIAss)), m_IETIPrec(gsScaledDirichletPrecondMPI<T>::make(*m_IETIAss))
        {
            m_IETISolv->init();
        }

        static uPtr make(typename gsIETIAssemblerMPI<T>::Ptr IETIAss) {return memory::make_unique(new IETIApplicationMPI(IETIAss));}

        virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
        {
            gsMatrix<T> temp,temp_acc;
            m_IETISolv->apply(input,temp);
            m_IETIPrec->accumulate(temp,temp_acc);
            m_IETIPrec->apply(temp_acc,x);
        }
        virtual void distribute(const gsMatrix<T>& input, gsMatrix<T> & distributed) const
        {
            m_IETIAss->distribute(input,distributed);
        }

        virtual void postAccumulate() const
        {
            m_IETIAss->postAccumulate();
        }
        virtual void startAccumulate(const  gsMatrix<T> & input) const
        {
            m_IETIAss->startAccumulate(input);
        }
        virtual void finishAccumulate(gsMatrix<T> & result) const
        {
            m_IETIAss->finishAccumulate(result);
        }

        /// Returns the number of rows of the operator
        virtual index_t rows() const {return m_IETISolv->rows();}

        /// Returns the number of columns of the operator
        virtual index_t cols() const {return rows();}

    private:
        typename gsIETIAssemblerMPI<T>::Ptr m_IETIAss;
        mutable typename gsIETISolverMPI<T>::Ptr m_IETISolv;
        mutable typename gsScaledDirichletPrecondMPI<T>::Ptr m_IETIPrec;
    };*/





    class ComplexSolverMPI : public gsLinearOperator<std::complex<T> >
    {
    public:

        /// Shared pointer for gsLinearOperator
        typedef memory::shared_ptr<ComplexSolverMPI> Ptr;

        /// Unique pointer for gsLinearOperator
        typedef memory::unique_ptr<ComplexSolverMPI> uPtr;

        ComplexSolverMPI(typename gsParallelOperator<T>::Ptr K, typename gsParallelOperator<T>::Ptr M, const std::complex<T>& lambda, typename gsDistributedOperator<T>::Ptr precond, gsMpiComm comm, int max_iter, real_t tol, bool prec, std::string outputFilename = "");

        static uPtr make(typename gsParallelOperator<T>::Ptr K, typename gsParallelOperator<T>::Ptr M, const std::complex<T>& lambda, typename gsDistributedOperator<T>::Ptr precond, gsMpiComm comm, int max_iter, real_t tol, bool prec, std::string outputFilename = "") {return memory::make_unique(new ComplexSolverMPI(K,M,lambda, precond,comm, max_iter, tol, prec,outputFilename));}

        virtual ~ComplexSolverMPI() {}

        /**
         * @brief apply the operator on the input vector and store the result in x
         * @param input Input vector
         * @param x     result vector
         */
        virtual void apply(const gsMatrix<std::complex<T> > & input, gsMatrix<std::complex<T> > & x) const;

        /// Returns the number of rows of the operator
        virtual index_t rows() const {return m_K->rows();}

        /// Returns the number of columns of the operator
        virtual index_t cols() const {return m_K->cols();}

    private:
        const typename gsParallelOperator<T>::Ptr m_K,m_M;
        std::complex<T> m_lambda;
        mutable gsMatrix<T > m_inp,  m_out;

        mutable gsMatrix<T > m_temp,  m_out_old;

        typename gsDistributedBlockOp<T>::Ptr blockOP;
        typename gsDistributedBlockOp<T>::Ptr blockPrec;

        int m_max_iter;
        real_t m_tol;
        bool m_print;
        std::string m_filename;


        typename gsParallelMinRes<T>::Ptr itSolver;
    };

    class RealSolverMPI : public gsLinearOperator<T >
    {
    public:

        /// Shared pointer for gsLinearOperator
        typedef memory::shared_ptr<RealSolverMPI> Ptr;

        /// Unique pointer for gsLinearOperator
        typedef memory::unique_ptr<RealSolverMPI> uPtr;


        RealSolverMPI(typename gsParallelOperator<T>::Ptr K, typename gsParallelOperator<T>::Ptr M, T a, T b, T c, typename gsDistributedOperator<T>::Ptr precond, int max_iter, real_t tol, std::string outputFilename = "");

        static uPtr make(typename gsParallelOperator<T>::Ptr K, typename gsParallelOperator<T>::Ptr M, T a, T b, T c, typename gsDistributedOperator<T>::Ptr precond, int max_iter, real_t tol, std::string outputFilename = "") {return memory::make_unique(new RealSolverMPI(K,M,a,b,c, precond, max_iter, tol,outputFilename));}

        virtual ~RealSolverMPI() {}

        /**
         * @brief apply the operator on the input vector and store the result in x
         * @param input Input vector
         * @param x     result vector
         */
        virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

        /// Returns the number of rows of the operator
        virtual index_t rows() const {return 2*m_K->rows();}

        /// Returns the number of columns of the operator
        virtual index_t cols() const {return 2*m_K->cols();}

    private:
        const typename gsParallelOperator<T>::Ptr m_K,m_M;
        T m_a,m_b,m_c;
        mutable gsMatrix<T > m_inp,  m_out;

        mutable gsMatrix<T > m_temp,  m_out_old;

        typename gsDistributedBlockOp<T>::Ptr blockOP;
        typename gsDistributedBlockOp<T>::Ptr blockPrec;

        int m_max_iter;
        real_t m_tol;
        bool m_print;
        std::string m_filename;

    //    typename gsGMRes<>::Ptr itSolver;
        typename gsParallelMinRes<real_t>::Ptr itSolver;
    };

protected:
    using Base::solverSetting;
    using Base::decompositionSetting;
    using Base::decompositionSolver;

    using Base::spaceLevels_;
    using Base::spaceLevel;
    using Base::ndofs_;

    using Base::NTimeLevels_;
    using Base::h_;

    using Base::transferSpace;
    using Base::transferTSpace;

    using Base::precA;
    using Base::LuOfA;

    using Base::settings_;

    using Base::tau_;
    using Base::MatTimeK_; using Base::MatTimeM_; using Base::MatTimeT0_;
    using Base::timeLevels_;

    using Base::myTimeSlices;

    using Base::spaceRestr_;
    using Base::timeRestr1_; using Base::timeRestr2_;

    using Base::A;
    using Base::B;
    using Base::data_;

    using Base::m_temp;
    mutable gsMatrix<T> m_temp_loc;

    using Base::TimeLevelsInp_;
    using Base::transferInp_;
    using Base::output_;

    gsSortedVector<size_t> myPatches_;
    gsMpiComm comm_;
    std::vector<typename gsPatchSubassembledTopology<T>::Ptr> m_subassTopology;
    std::vector<gsParallelGlobalLocalHandler::Ptr> m_parSpaceHandler;
    gsParallelGlobalLocalHandler::Ptr m_parHandler;
    std::vector<typename gsParallelOperator<T>::Ptr > m_prolongation;
    std::vector<typename gsParallelOperator<T>::Ptr > m_restriction;

    std::vector<typename gsParallelOperator<T>::Ptr> m_Kx;
    std::vector<typename gsParallelOperator<T>::Ptr> m_Mx;
    std::vector<typename gsParallelOperator<T>::Ptr> m_MWx;

    std::vector<typename gsMatrixOp< gsSparseMatrix<T> >::Ptr > m_coarseKx;
    std::vector<typename gsMatrixOp< gsSparseMatrix<T> >::Ptr > m_coarseMx;
    std::vector<gsDofMapper> m_coarseLocMappers;
    gsDofMapper m_coarseGlobMapper;

public:
    gsSpaceParallelHSTSlapOp(real_t h, real_t tau, int NTimeLevels);

    /// Make function returning a smart pointer
    static uPtr make(real_t h, real_t tau, int NTimeLevels){ return memory::make_unique( new gsSpaceParallelHSTSlapOp(h, tau,NTimeLevels) ); }

    //Matrices are already in Tensor form
    gsSpaceParallelHSTSlapOp(const std::vector<int> &TimeLevels, int NTimeLevels,
                             int spaceLevels,
                             int spaceLevel,
                             std::vector<std::vector<int> > myTimeSlices,
                             const gsSortedVector<size_t>&  myPatches,
                             std::vector<STData > data,
                             std::vector<typename gsParallelOperator<T>::Ptr > prolongation,
                             std::vector<typename gsParallelOperator<T>::Ptr > restriction,
                             std::vector<typename gsPatchSubassembledTopology<T>::Ptr> subassTopology,
                             std::vector<typename gsParallelGlobalLocalHandler::Ptr> handlers,
                             std::vector<typename gsParallelOperator<T>::Ptr> Kx,
                             std::vector<typename gsParallelOperator<T>::Ptr> Mx,
                             std::vector<typename gsParallelOperator<T>::Ptr> Mwx,
                             gsMatrix<T> timeRestr1,
                             gsMatrix<T> timeRestr2,
                             real_t tau,
                             real_t h,
                             typename Base::HeatSTSlapSettings settings, const gsMpiComm& comm, bool output=true);

    static uPtr make(const std::vector<int> &TimeLevels, int NTimeLevels,
                     int spaceLevels,
                     int spaceLevel,
                     std::vector<std::vector<int> > myTimeSlices,
                     const gsSortedVector<size_t>&  myPatches,
                     std::vector<STData> data,
                     std::vector<typename gsParallelOperator<T>::Ptr > prolongation,
                     std::vector<typename gsParallelOperator<T>::Ptr > restriction,
                     std::vector<typename gsPatchSubassembledTopology<T>::Ptr> subassTopology,
                     std::vector<typename gsParallelGlobalLocalHandler::Ptr> handlers,
                     std::vector<typename gsParallelOperator<T>::Ptr> Kx,
                     std::vector<typename gsParallelOperator<T>::Ptr> Mx,
                     std::vector<typename gsParallelOperator<T>::Ptr> Mwx,
                     gsMatrix<T> timeRestr1,
                     gsMatrix<T> timeRestr2,
                     real_t tau,
                     real_t h,
                     typename Base::HeatSTSlapSettings settings,const gsMpiComm& comm, bool output=true)
    { return memory::make_unique( new gsSpaceParallelHSTSlapOp(TimeLevels, NTimeLevels,spaceLevels,spaceLevel,myTimeSlices,myPatches,data,prolongation,restriction,subassTopology,handlers,Kx,Mx,Mwx,timeRestr1,timeRestr2,tau,h,settings,comm,output) ); }

    virtual void init();

    virtual void setMGCoarseMatrices(std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr > coarseKx,
                                     std::vector<typename gsMatrixOp<gsSparseMatrix<T> >::Ptr > coarseMx,
                                     std::vector<gsDofMapper> coarseLocMappers,
                                     gsDofMapper coarseGlobMapper);

    virtual void calcInitialVector(const Vector &u, Vector &f, int timeLevel, int timeStep) const;

    virtual void mult(const Vector &u, Vector &f, int timeLevel, int timeStep, real_t sign) const;

    void extractAndReorderRhs(gsMatrix<T>& rhs) const;
    void buildSliceWiseRHS(std::vector<gsMatrix<T> >& rhs, const std::vector<int> & sliceIdx)const;
    void buildSliceWiseSolution(std::vector< gsMatrix<T> > & sol, const std::vector<int> & sliceIdx) const;
    void accumulate(const Vector & v,Vector & v_acc, MPI_Comm spaceComm) const;
    void distribute(const Vector & v,Vector & v_acc, MPI_Comm spaceComm) const;

    virtual void TimeRestriction1(const Vector &u_fine, Vector &u_coarse, int coarseTimeLevel, bool setZero=true) const;
    virtual void TimeRestriction2(const Vector &u_fine, Vector &u_coarse, int coarseTimeLevel, bool setZero=true) const;

    virtual void TimeProlongation1(const Vector &u_coarse, Vector &u_fine, int coarseTimeLevel, bool setZero=true) const;
    virtual void TimeProlongation2(const Vector &u_coarse, Vector &u_fine, int coarseTimeLevel, bool setZero=true) const;


};

} // namespace gismo
#endif

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMGSpaceTimeAdapter_MPI.hpp)
#endif

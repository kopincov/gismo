/** @file gsMGSpaceTimeAdapter.h

    @brief Adapter Class for gsTimeParallelMultigrid class  considering the heat equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
    Created on: 2017-11-23
*/


#pragma once
#include <gsSolver/gsMinimalResidual.h>
#include <gsSolver/gsBlockOp.h>
#include <gsAssembler/gsSpaceTimesliceAssembler.h>


namespace gismo {


template<typename T> class gsParallelOperator;

namespace Smoother_ {
enum type {
    Richardson,
    Jacobi,
    GaussSeidel,
    MassRichardsonSubspaceCorrectionAdditiveMPDDD,
    MassRichardsonSubspaceCorrectionAdditiveMPDDD2,
};
}

template<typename T>
class gsHeatSTSlapOperator
{

public:
    //#define localMatrix BlasMatrix
    typedef  gsMatrix<T> localMatrix;
    typedef gsMatrix<T> Vector;

    /// Shared pointer for gsIdentityOp
    typedef memory::shared_ptr<gsHeatSTSlapOperator> Ptr;

    /// Unique pointer for gsIdentityOp
    typedef memory::unique_ptr<gsHeatSTSlapOperator> uPtr;

    //typedef typename gsSpaceTimesliceAssembler<T>::TPData TPData;
    struct STData : public gsSpaceTimesliceAssembler<T>::TPData
    {
        std::vector<memory::shared_ptr<gsParallelOperator<T> > > parKx, parMx, parMWx;
    };

    enum solverSetting
    {
        FULLSPACETIME = 0,
        SPATIALWISE,
    };

    enum decompositionSolver
    {
        DIRECT=0,
        IETI,
        MG
    };

    enum decompositionSetting
    {
        DIAGONALIZATION=0,
        REALSCHUR,
        COMPLEXSCHUR
    };

    struct HeatSTSlapSettings
    {
        HeatSTSlapSettings(gsOptionList opt);

        HeatSTSlapSettings();

        static gsOptionList defaultOptions();

        Smoother_::type smoother;
        real_t damping;
        real_t outerDamping;
        int sliceSolver;
        int spaceSolver;
        int timeDecomp;
        int approxIterates;
        int precondIterations;
        real_t precondTol;
        real_t tol;
        gsOptionList IETIOptions;
        gsOptionList MGOptions;
        std::vector<memory::shared_ptr<gsAssembler<T> > > m_spatialAssemblers; //required for IETI and MG;
        std::vector<memory::shared_ptr<gsAssembler<T> > > m_STAssemblers; //required for IETI;
        std::vector<memory::shared_ptr<gsPde<T> > > m_pde;
    };


    class ComplexSolver : public gsLinearOperator<std::complex<T> >
    {
    public:

        /// Shared pointer for gsLinearOperator
        typedef memory::shared_ptr<ComplexSolver> Ptr;

        /// Unique pointer for gsLinearOperator
        typedef memory::unique_ptr<ComplexSolver> uPtr;

        ComplexSolver(const typename gsSparseMatrix<T>::Ptr& K, const  typename gsSparseMatrix<T>::Ptr& M, const std::complex<T>& lambda, typename gsLinearOperator<T>::Ptr precond, int max_iter, real_t tol, bool prec, std::string outputFilename = "");

        static uPtr make(const typename gsSparseMatrix<T>::Ptr& K, const  typename gsSparseMatrix<T>::Ptr& M, const std::complex<T>& lambda, typename gsLinearOperator<T>::Ptr precond, int max_iter, real_t tol, bool prec, std::string outputFilename = "") {return memory::make_unique(new ComplexSolver(K,M,lambda, precond, max_iter, tol, prec,outputFilename));}

        virtual ~ComplexSolver() {}

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
        const typename gsSparseMatrix<T>::Ptr m_K,m_M;
        std::complex<T> m_lambda;
        mutable gsMatrix<T > m_inp,  m_out;

        mutable gsMatrix<T > m_temp,  m_out_old;

        typename gsBlockOp<T>::Ptr blockOP;
        typename gsBlockOp<T>::Ptr blockPrec;

        int m_max_iter;
        real_t m_tol;
        bool m_print;
        std::string m_filename;


        typename gsMinimalResidual<>::Ptr itSolver;
    };

    class RealSolver : public gsLinearOperator<T >
    {
    public:

        /// Shared pointer for gsLinearOperator
        typedef memory::shared_ptr<RealSolver> Ptr;

        /// Unique pointer for gsLinearOperator
        typedef memory::unique_ptr<RealSolver> uPtr;


        RealSolver(const typename gsSparseMatrix<T>::Ptr& K, const  typename gsSparseMatrix<T>::Ptr& M, T a, T b, T c, typename gsLinearOperator<T>::Ptr precond, int max_iter, real_t tol, std::string outputFilename = "");

        static uPtr make(const typename gsSparseMatrix<T>::Ptr& K, const  typename gsSparseMatrix<T>::Ptr& M, T a, T b, T c, typename gsLinearOperator<T>::Ptr precond, int max_iter, real_t tol, std::string outputFilename = "") {return memory::make_unique(new RealSolver(K,M,a,b,c, precond, max_iter, tol,outputFilename));}

        virtual ~RealSolver() {}

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
        const typename gsSparseMatrix<T>::Ptr m_K,m_M;
        T m_a,m_b,m_c;
        mutable gsMatrix<T > m_inp,  m_out;

        mutable gsMatrix<T > m_temp,  m_out_old;

        typename gsBlockOp<T>::Ptr blockOP;
        typename gsBlockOp<T>::Ptr blockPrec;

        int m_max_iter;
        real_t m_tol;
	bool m_print;
	std::string m_filename;

    //    typename gsGMRes<>::Ptr itSolver;
          typename gsMinimalResidual<real_t>::Ptr itSolver;
    };


    //optimized class for the special matrix structure:
    //  [A1 tB tB tB]
    //  [   A2 tB tB]
    //  [      .. tB]
    //  [         AN]
    //The scalar value t may be different in each block, however, B is the same.


    class BlockTriangularSolver : public gsLinearOperator<std::complex<T> >
    {
    public:

        /// Shared pointer for gsLinearOperator
        typedef memory::shared_ptr<BlockTriangularSolver> Ptr;

        /// Unique pointer for gsLinearOperator
        typedef memory::unique_ptr<BlockTriangularSolver> uPtr;

        BlockTriangularSolver(std::vector<typename gsLinearOperator<std::complex<T> >::Ptr> A_inv, typename gsSparseMatrix<T>::Ptr B, typename gsMatrix<std::complex<T> >::Ptr UT);
        BlockTriangularSolver(std::vector<typename gsLinearOperator<std::complex<T> >::Ptr> A_inv, typename gsLinearOperator<T>::Ptr B, typename gsMatrix<std::complex<T> >::Ptr UT);


        static uPtr make(std::vector<typename gsLinearOperator<std::complex<T> >::Ptr> A_inv, typename gsSparseMatrix<T>::Ptr B, typename gsMatrix<std::complex<T> >::Ptr UT) {return memory::make_unique(new BlockTriangularSolver(A_inv,B,UT));}
        static uPtr make(std::vector<typename gsLinearOperator<std::complex<T> >::Ptr> A_inv, typename gsLinearOperator<T>::Ptr B, typename gsMatrix<std::complex<T> >::Ptr UT) {return memory::make_unique(new BlockTriangularSolver(A_inv,B,UT));}


        virtual ~BlockTriangularSolver() {}

        /**
         * @brief apply the operator on the input vector and store the result in x
         * @param input Input vector
         * @param x     result vector
         */
        virtual void apply(const gsMatrix<std::complex<T> > & input, gsMatrix<std::complex<T> > & x) const;

        /// Returns the number of rows of the operator
        virtual index_t rows() const {return m_Ainv.size()*m_B->rows();}

        /// Returns the number of columns of the operator
        virtual index_t cols() const {return rows();}

    private:
        std::vector<typename gsLinearOperator<std::complex<T> >::Ptr> m_Ainv;
        typename gsLinearOperator<T>::Ptr m_B;
        typename gsMatrix<std::complex<T> >::Ptr m_T;

        mutable gsVector<index_t> m_rowVector, m_colVector;
        mutable std::vector<gsMatrix<std::complex<T> > > m_z;
        mutable gsMatrix<std::complex<T> > m_temp,  m_rhs;
        mutable gsMatrix<T> m_realTemp;

    };

    class QuasiBlockTriangularSolver : public gsLinearOperator<T >
    {
    public:

        /// Shared pointer for gsLinearOperator
        typedef memory::shared_ptr<QuasiBlockTriangularSolver> Ptr;

        /// Unique pointer for gsLinearOperator
        typedef memory::unique_ptr<QuasiBlockTriangularSolver> uPtr;

        QuasiBlockTriangularSolver(std::vector<typename gsLinearOperator<T>::Ptr> A_inv, typename gsSparseMatrix<T>::Ptr B, typename gsMatrix<T>::Ptr UT);
        QuasiBlockTriangularSolver(std::vector<typename gsLinearOperator<T>::Ptr> A_inv, typename gsLinearOperator<T>::Ptr B, typename gsMatrix<T>::Ptr UT);

        static uPtr make(std::vector<typename gsLinearOperator<T>::Ptr> A_inv, typename gsSparseMatrix<T>::Ptr B, typename gsMatrix<T>::Ptr UT) {return memory::make_unique(new QuasiBlockTriangularSolver(A_inv,B,UT));}
        static uPtr make(std::vector<typename gsLinearOperator<T>::Ptr> A_inv, typename gsLinearOperator<T>::Ptr B, typename gsMatrix<T>::Ptr UT) {return memory::make_unique(new QuasiBlockTriangularSolver(A_inv,B,UT));}

        virtual ~QuasiBlockTriangularSolver() {}

        /**
         * @brief apply the operator on the input vector and store the result in x
         * @param input Input vector
         * @param x     result vector
         */
        virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

        /// Returns the number of rows of the operator
        virtual index_t rows() const {return m_T->rows()*m_B->rows();}

        /// Returns the number of columns of the operator
        virtual index_t cols() const {return rows();}

    private:
        std::vector<typename gsLinearOperator<T>::Ptr> m_Ainv;
        typename gsLinearOperator<T>::Ptr m_B;
        typename gsMatrix<T>::Ptr m_T;

        mutable gsVector<index_t> m_rowVector, m_colVector;
        mutable std::vector<gsMatrix<T> > m_z;
        mutable gsMatrix<T> m_temp,  m_rhs;

    };

    /*
    class IETIApplication : public gsLinearOperator<T>
    {
    public:
        /// Shared pointer for gsLinearOperator
        typedef memory::shared_ptr<IETIApplication> Ptr;

        /// Unique pointer for gsLinearOperator
        typedef memory::unique_ptr<IETIApplication> uPtr;


        IETIApplication(typename gsIETIAssembler<T>::Ptr IETIAss): m_IETIAss(IETIAss) , m_IETISolv(gsIETISolver<T>::make(*m_IETIAss)), m_IETIPrec(gsScaledDirichletPrecond<T>::make(*m_IETIAss))
        {
            m_IETISolv->init();
        }

        static uPtr make(typename gsIETIAssembler<T>::Ptr IETIAss) {return memory::make_unique(new IETIApplication(IETIAss));}

        virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
        {
            gsMatrix<T> temp;
            m_IETISolv->apply(input,temp);
            m_IETIPrec->apply(temp,x);
        }

        /// Returns the number of rows of the operator
        virtual index_t rows() const {return m_IETISolv->rows();}

        /// Returns the number of columns of the operator
        virtual index_t cols() const {return rows();}

    private:
        typename gsIETIAssembler<T>::Ptr m_IETIAss;
        mutable typename gsIETISolver<T>::Ptr m_IETISolv;
        mutable typename gsScaledDirichletPrecond<T>::Ptr m_IETIPrec;
    };
*/

/*
    template<bool inverse>
    class gsApplyPermutation : public gsLinearOperator<T>
    {
    public:

        /// Shared pointer for gsLinearOperator
        typedef memory::shared_ptr<gsApplyPermutation> Ptr;

        /// Unique pointer for gsLinearOperator
        typedef memory::unique_ptr<gsApplyPermutation> uPtr;

        gsApplyPermutation(memory::shared_ptr<Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic, index_t> > perm) :m_perm(perm) {}
        static uPtr make(memory::shared_ptr<Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic, index_t> > perm) {return memory::make_unique(new gsApplyPermutation(perm));}

        virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
        {
            if(inverse)
                x = m_perm->inverse()*input;
            else
                x = (*m_perm)*input;
        }

        /// Returns the number of rows of the operator
        virtual index_t rows() const {return m_perm->rows();}

        /// Returns the number of columns of the operator
        virtual index_t cols() const {return rows();}

    private:
        memory::shared_ptr<Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic, index_t> > m_perm;
    };
*/

private:
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


protected:
    //  BlockMatrix<localMatrix> M0_; SparseMatrix tempM0_;
    // SparseMatrix *Mh_;
    int spaceLevels_;
    int spaceLevel;
    // int m_, n_;
    //int ndofs_;
    std::vector<int> ndofs_;

    int NTimeLevels_;
    real_t h_;

    // std::vector<std::vector<BlockMatrix<localMatrix> > > A_;
    //std::vector<typename gsKroneckerOp<T>::uPtr> transferTime;
    std::vector<typename gsKroneckerOp<T>::Ptr> transferSpace;
    std::vector<typename gsKroneckerOp<T>::Ptr> transferTSpace;
    //  std::vector<BlockMG<CoarseSolver, localMatrix>*> GMG_;

    std::vector<std::vector<typename  gsLinearOperator<T>::Ptr > > precA;
    std::vector<std::vector<typename  gsLinearOperator<T>::Ptr > > LuOfA;

    //multigrid settings
    HeatSTSlapSettings settings_;

    real_t tau_;
    //  int order_;
    //  int blocksize_;
    std::vector<gsMatrix<T> > MatTimeK_, MatTimeM_, MatTimeT0_;
    std::vector<int> timeLevels_;

    std::vector<std::vector<int> > myTimeSlices;

    std::vector<typename gsSparseMatrix<T,RowMajor>::Ptr > spaceRestr_;
    gsMatrix<T> timeRestr1_, timeRestr2_;

    std::vector<std::vector<typename gsLinearOperator<T>::Ptr > > A;
    std::vector<std::vector<typename gsLinearOperator<T>::Ptr > > B;

    // std::vector<gsLinearOperator::uPtr > A;
    // std::vector<gsLinearOperator::uPtr > B;

    std::vector<STData > data_;

    mutable gsMatrix<T> m_temp;

    std::vector<int> TimeLevelsInp_;
    std::vector<typename gsSparseMatrix<T,RowMajor>::Ptr > transferInp_;
    bool output_;
    std::string m_filename;

public:
    gsHeatSTSlapOperator(real_t h, real_t tau, int NTimeLevels);

    /// Make function returning a smart pointer
    static uPtr make(real_t h, real_t tau, int NTimeLevels){ return memory::make_unique( new gsHeatSTSlapOperator(h, tau,NTimeLevels) ); }
    /// Make function returning a smart pointer

    //gsMatrix<real_t> getMassTMatrix(int timeLevel) const { return MatTimeM_[timeLevels_[timeLevel]]; }

    //Matrices are already in Tensor form
    gsHeatSTSlapOperator(const std::vector<int> &TimeLevels, int NTimeLevels,
                         int spaceLevels,
                         int spaceLevel,
                         std::vector<std::vector<int> > myTimeSlices,
                         std::vector<STData > data,
                         std::vector<typename gsSparseMatrix<T,RowMajor>::Ptr >& transfer,
                         gsMatrix<T> timeRestr1,
                         gsMatrix<T> timeRestr2,
                         real_t tau,
                         real_t h,
                         HeatSTSlapSettings settings, bool output=true);

    static uPtr make(const std::vector<int> &TimeLevels, int NTimeLevels,
                     int spaceLevels,
                     int spaceLevel,
                     std::vector<std::vector<int> > myTimeSlices,
                     std::vector<STData> data,
                     std::vector<typename gsSparseMatrix<T,RowMajor>::Ptr >& transfer,
                     gsMatrix<T> timeRestr1,
                     gsMatrix<T> timeRestr2,
                     real_t tau,
                     real_t h,
                     HeatSTSlapSettings settings, bool output=true)
    { return memory::make_unique( new gsHeatSTSlapOperator(TimeLevels, NTimeLevels,spaceLevels,spaceLevel,myTimeSlices,data,transfer,timeRestr1,timeRestr2,tau,h,settings,output) ); }
protected:
    gsHeatSTSlapOperator(const std::vector<int> &TimeLevels, int NTimeLevels,
                         int spaceLevels,
                         int spaceLevel,
                         std::vector<std::vector<int> > myTimeSlices,
                         std::vector<STData > data,
                         gsMatrix<T> timeRestr1,
                         gsMatrix<T> timeRestr2,
                         real_t tau,
                         real_t h,
                         HeatSTSlapSettings settings, bool output=true);


public:


    virtual void init();

    void setOutputFilename(std::string filename) { m_filename = filename; }

    virtual int solveOnSTSlap(Vector &u, const Vector &f, int timeLevel, int timeStep) const;

    virtual void calcInitialVector(const Vector &u, Vector &f, int timeLevel, int timeStep) const;

    //TODO: identicalSlice
    int getNDofs(void) const ;

    //TODO: identicalSlice
    bool checkIfSpaceCoarseningIsAllowed(int timeLevel) const;

    virtual void mult(const Vector &u, Vector &f, int timeLevel, int timeStep, real_t sign) const;

    //TODO: identicalSlice
    virtual void ApproximateSolve(Vector &u, const Vector &f, int timeLevel, int timeStep) const;

    //Restrictions are time-slice wise
    //Two timeslices are restricted to one
    //On a timeslice the basis does not change (if all timeslices have the same basis).


    //TODO: identicalSlice
    virtual void TimeRestriction1(const Vector &u_fine, Vector &u_coarse, int coarseTimeLevel, bool setZero=true) const;
    //TODO: identicalSlice
    virtual void TimeRestriction2(const Vector &u_fine, Vector &u_coarse, int coarseTimeLevel, bool setZero=true) const;
    //TODO: identicalSlice
    virtual void SpaceRestriction(const Vector &u_fine, Vector &u_coarse, int coarseSpaceLevel, bool setZero=true) const;
    //TODO: identicalSlice
    virtual void TimeProlongation1(const Vector &u_coarse, Vector &u_fine, int coarseTimeLevel, bool setZero=true) const;

    //TODO:identicalSlice
    virtual void TimeProlongation2(const Vector &u_coarse, Vector &u_fine, int coarseTimeLevel, bool setZero=true) const;

    //TODO: identicalSlice
    virtual void SpaceProlongation(const Vector &u_coarse, Vector &u_fine, int coarseSpaceLevel, bool setZero=true) const;
};








} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMGSpaceTimeAdapter.hpp)
#endif

/** @file gsBramblePasciakCG.h

    @brief General Conjugate gradient solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer
*/

#pragma once

#include <gsSolver/gsIterativeSolver.h>
#include <gsSolver/gsBlockOp.h>

#include <gsSolver/gsSolverUtils.h>

namespace gismo
{

namespace gsBPCG_Types
{
/// Type for \gsBramblePasciakCG 
enum Type
{
    BP = 0,                ///< Stanadard
    BP_Schur = 1,          ///< Schur
    SchoeberlZulehner = 2, ///< SchoeberlZulehner
};
} // namespace gsBPCG_Types

template<gsBPCG_Types::Type type, typename T> class gsBramblePasciakCG;

namespace internal
{

template<gsBPCG_Types::Type type, typename VectorType> struct gsBPCG_Preconditioner{};

template<typename VectorType>
struct gsBPCG_Preconditioner<gsBPCG_Types::BP,VectorType>
{
    typedef gsVector<typename VectorType::Scalar,1> CombinationScalingType;
    typedef gsVector<typename VectorType::Scalar,1> PreconditionScalingType;
    static void apply( VectorType& input, VectorType& result, VectorType& hResult,
                       gsBramblePasciakCG<gsBPCG_Types::BP,typename VectorType::Scalar>& bpcg );
};

template<typename VectorType>
struct gsBPCG_Preconditioner<gsBPCG_Types::BP_Schur,VectorType>
{
    typedef gsVector<typename VectorType::Scalar,1> CombinationScalingType;
    typedef gsVector<typename VectorType::Scalar,2> PreconditionScalingType;
    static void apply( VectorType& input, VectorType& result, VectorType& hResult,
                       gsBramblePasciakCG<gsBPCG_Types::BP_Schur,typename VectorType::Scalar>& bpcg );
};

template<typename VectorType>
struct gsBPCG_Preconditioner<gsBPCG_Types::SchoeberlZulehner,VectorType>
{
    typedef gsVector<typename VectorType::Scalar,0> CombinationScalingType;
    typedef gsVector<typename VectorType::Scalar,2> PreconditionScalingType;
    static void apply( VectorType& input, VectorType& result, VectorType& hResult,
                       gsBramblePasciakCG<gsBPCG_Types::SchoeberlZulehner,typename VectorType::Scalar>& bpcg );
};

} // namespace internal

template<gsBPCG_Types::Type type = gsBPCG_Types::BP, typename T=real_t>
class gsBramblePasciakCG : public gsIterativeSolver<T>
{
protected:
    
    template <gsBPCG_Types::Type,typename> friend struct internal::gsBPCG_Preconditioner;

public:

    typedef gsIterativeSolver<T> Base;

    typedef gsMatrix<T> VectorType;

    typedef typename internal::gsBPCG_Preconditioner<type,VectorType>::CombinationScalingType CombinationScalingType;

    typedef typename internal::gsBPCG_Preconditioner<type,VectorType>::PreconditionScalingType PreconditionScalingType;
    
    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsBramblePasciakCG<type,T> > Ptr;
    typedef memory::unique_ptr<gsBramblePasciakCG<type,T> > uPtr;

    /// @brief Constructor using a matrix (operator) and optionally preconditionners
    ///
    /// @param mat                     The operator to be solved for; here we require a gsBlockOp
    /// @param precondA                The preconditioner for the (1,1)-block, defaulted to the identity
    /// @param precondC                The preconditioner for the (2,2)-block, defaulted to the identity
    /// @param CombinationScalingType  A scaling type constant
    /// @param PreconditionScalingType A scaling type constant
    
    explicit gsBramblePasciakCG( const typename gsBlockOp<T>::Ptr& mat, const LinOpPtr & precondA = LinOpPtr(), const LinOpPtr & precondC = LinOpPtr(),
                                 const CombinationScalingType& scaling = CombinationScalingType::Constant(1),
                                 const PreconditionScalingType& precondScaling = PreconditionScalingType::Constant(1) )
        : Base(mat, LinOpPtr()), m_precondA(precondA), m_precondC(precondC), m_scaling(scaling), m_precond_scaling(precondScaling), m_calcEigenvals(false)
    {
        
        GISMO_ASSERT( m_precondA,
                      "At least the preconditioner for A must be defined." );

        GISMO_ENSURE( mat->rowBlocks()==2 && mat->colBlocks()==2,
                      "Method can handle only 2x2 block systems." );
        
        m_matA = mat->getOperator(0,0);
        m_matB = mat->getOperator(1,0);
        m_matBT = mat->getOperator(0,1);
        //m_matC = blOp->getOperator(1,1);

        if (!m_precondC) m_precondC = gsIdentityOp<T>::make(m_mat->rows()-m_precondA->rows());

        // Consistency of the matrices is guaranteed by gsBlockOp. We just need to ensure that A and C are square.
        // Note that C needs not to be defined.

        GISMO_ASSERT( m_matA->rows() == m_matA->cols() && m_matB->rows() == m_matBT->cols(),
                      "The matrices A and C need to be square.");
        
        GISMO_ASSERT( m_precondA->rows() == m_precondA->cols() && m_precondA->rows() == m_matA->rows(),
                      "The preconditioner for A is not square or does not match the matrix." );

        GISMO_ASSERT( m_precondC->rows() == m_precondC->cols() && m_precondC->rows() == m_matB->rows(),
                      "The preconditioner for C is not square or does not match the matrix." );
        
        m_nA = m_matA->rows();
        m_nC = m_matB->rows();
    }

    /// @brief Make function using a matrix (operator) and optionally preconditionners
    ///
    /// @param mat                     The operator to be solved for; here we require a gsBlockOp
    /// @param precondA                The preconditioner for the (1,1)-block, defaulted to the identity
    /// @param precondC                The preconditioner for the (2,2)-block, defaulted to the identity
    /// @param CombinationScalingType  A scaling type constant
    /// @param PreconditionScalingType A scaling type constant
    static uPtr make( const typename gsBlockOp<T>::Ptr& mat, const LinOpPtr & precondA = LinOpPtr(), const LinOpPtr & precondC = LinOpPtr(),
                      const CombinationScalingType& scaling = CombinationScalingType::Constant(1),
                      const PreconditionScalingType& precondScaling = PreconditionScalingType::Constant(1) )
    { return uPtr( new gsBramblePasciakCG<type,T>( mat, precondA, precondC,scaling,precondScaling) ); }
    
    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addSwitch("CalcEigenvalues", "Additionally to solving the system,"
                                         " CG computes the eigenvalues of the Lanczos matrix", false );
        if(CombinationScalingType::RowsAtCompileTime>0)
            opt.addReal("CombinationScaling", "Implemented for BP and BP_Schur, may improve performance", 1);
        if(PreconditionScalingType::RowsAtCompileTime>0)
            opt.addReal("PreconditionScaling1", "First parameter for scaling the preconditioner, required for matrix inequality",1);
        if(PreconditionScalingType::RowsAtCompileTime>1)
            opt.addReal("PreconditionScaling2", "Second scaling parameter (required by SchöberlZulehner)",1);
        return opt;
    }

    /// @brief Set the options based on a gsOptionList
    gsBramblePasciakCG& setOptions( const gsOptionList & opt )
    {
        Base::setOptions(opt);
        m_calcEigenvals = opt.askSwitch("CalcEigenvalues", m_calcEigenvals);
        if(m_scaling.RowsAtCompileTime>0)
            m_scaling(0) = opt.askReal("CombinationScaling", m_scaling(0));
        if(m_precond_scaling.RowsAtCompileTime>0)
            m_precond_scaling(0) = opt.askReal("PreconditionScaling1", m_precond_scaling(0));
        if(m_precond_scaling.RowsAtCompileTime>1)
            m_precond_scaling(1) = opt.askReal("PreconditionScaling2", m_precond_scaling(1));

        return *this;
    }

    void applyPreconditioner( VectorType& input, VectorType& result, VectorType& hResult )
    {
        internal::gsBPCG_Preconditioner<type,VectorType>::apply(input,result,hResult,*this);
    }

    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );


    /// @brief Specify if you want to store data for eigenvalue estimation
    /// @param flag true stores the coefficients of the lancos matrix, false not.
    void setCalcEigenvalues( bool flag )       { m_calcEigenvals = flag;      }
    
    
    /// @brief Specify the scaling parameter TODO: which one?
    void setCombinationScaling(const CombinationScalingType & scaling )               { m_scaling = scaling;                                              }
    void setCombinationScaling(T scaling )                                            { m_scaling = CombinationScalingType::Constant(scaling);            }

    /// @brief Specify the scaling parameter TODO: which one?
    void setPreconditionerScaling(const PreconditionScalingType& scaling ) { m_precond_scaling = scaling;                                      }
    void setPreconditionerScaling(T scaling )                              { m_precond_scaling = PreconditionScalingType::Constant(scaling);   }

    /// @brief returns the condition number of the (preconditioned) system matrix
    T getConditionNumber();

    /// @brief returns the eigenvalues of the Lanczos matrix
    void getEigenvalues( gsMatrix<T>& eigs );

    /// Prints the object as a string.
    std::ostream &print( std::ostream &os ) const
    {
        os << "gsBramblePasciakCG\n";
        return os;
    }

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;

    LinOpPtr m_matA;
    LinOpPtr m_matB;
    LinOpPtr m_matBT;
    //LinOpPtr m_matC;

    LinOpPtr m_precondA;
    LinOpPtr m_precondC;

    CombinationScalingType   m_scaling;
    PreconditionScalingType  m_precond_scaling;

    VectorType m_z;
    VectorType m_hz;
    VectorType m_update;
    VectorType m_tmp;
    VectorType m_hPtmp;
    VectorType m_precTemp;
    VectorType m_residual;
    T m_abs_new;


    bool m_calcEigenvals;

    index_t m_nA;
    index_t m_nC;

    std::vector<T> m_delta, m_gamma;
};

} // namespace gismo


//#ifndef GISMO_BUILD_LIB
//#include GISMO_HPP_HEADER(gsBramblePasciakCG.hpp)
//#endif

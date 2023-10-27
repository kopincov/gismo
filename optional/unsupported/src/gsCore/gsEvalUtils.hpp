/** @file gsEvalUtils.hpp

    @brief Functions that compute div, jacobian, and curl from g+smo deriv format
    and hessian and laplacian from the deriv2 format.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, A. Mantzaflaris
**/

#pragma once

// TODO test, clean up, optimize
// one possibility is to use Eigen::Map of proper size and orientation to avoid copying
// maybe it is worth to template the methods on the domain and target dimension and dispetch
// to the specialized template on the base class (leaving the common implementation as default)

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunctionSet.h>

namespace gismo {



template <typename T>
struct convertValue
{
    static inline  void derivToJacobian  (const gsMatrix<T> &deriv,
                                          const std::pair<short_t, short_t> & Dim,
                                          gsMatrix<T> &result)
    {
        const int blockSize = Dim.first * Dim.second;
        const int numBlocks = deriv.rows() / blockSize; // num. functions

        result.resize(numBlocks*Dim.second,Dim.first*deriv.cols());
        for (int p=0; p<deriv.cols(); ++p)
            for (int b=0; b< numBlocks; ++b)
            {
                gsMatrix<T> temp2=deriv.block(blockSize*b,p,blockSize,1);
                temp2.resize(Dim.second,Dim.first);
                result.block(b*Dim.second,p*Dim.first,Dim.second,Dim.first)=temp2;
            }
    }

    static inline  void derivToDiv  (const gsMatrix<T> &deriv,
                                     const std::pair<short_t, short_t> & Dim,
                                     gsMatrix<T> &result)
    {
        GISMO_ASSERT( (Dim.second%Dim.first)==0,
                      "Incompatible domain and target for divergence");
        const int blockSize = Dim.first * Dim.second;
        const int numBlocks = deriv.rows() / blockSize; // num. functions
        //gsInfo << "deriv : " << deriv << "\n";
        result.setZero(numBlocks,deriv.cols());
        for (int b=0; b< numBlocks; ++b)
            for (int s=0; s<Dim.first;++s)
                result.row(b)+=deriv.row(b*blockSize+s*(Dim.first+1));
        //gsInfo << "derivToDiv : " << result << "\n";
    }

    static inline  void derivToCurlOne(
            typename gsMatrix<T>::constRows deriv,
            typename gsMatrix<T>::Rows      curl,
            short_t domDim
            )
    {
        switch ( domDim )
        {
        case 1:
            curl=deriv;
            break;
        case 2:
            curl.row(0)=-deriv.row(1);
            curl.row(1)=deriv.row(0);
            break;
        case 3:
            curl.row(0)=deriv.row(7)-deriv.row(5);
            curl.row(1)=deriv.row(2)-deriv.row(6);
            curl.row(2)=deriv.row(3)-deriv.row(1);
            break;
        default:
            GISMO_ERROR("curl is not implemented in 4d  or more");
        }
    }

    static inline  void derivToCurl  (const gsMatrix<T> &deriv,
                                      const std::pair<short_t, short_t> & Dim,
                                      gsMatrix<T> &result)
    {
        int outBlockSize;
        int inBlockSize;
        switch (Dim.first)
        {
        case 1:
            inBlockSize  = 1;
            outBlockSize = 1;
            break;
        case 2:
            inBlockSize  = 2;
            outBlockSize = 2;
            break;
        case 3:
            inBlockSize  = 9;
            outBlockSize = 3;
            break;
        default:
            GISMO_ERROR("no implementation of curl for this domain dimension");
        }
        GISMO_ASSERT( (Dim.first*Dim.second)%inBlockSize == 0, "incompatible domain and target dimension for curl operator" );

        const int numBlocks = deriv.rows()/inBlockSize;
        const short_t domDim    = Dim.first;

        result.resize(numBlocks*outBlockSize,deriv.cols());
        for (int b=0; b< numBlocks; ++b)
            derivToCurlOne(
                        deriv.middleRows(inBlockSize*b,inBlockSize),
                        result.middleRows(outBlockSize*b,outBlockSize),
                        domDim);
    }


    static inline  void deriv2ToLaplacian  (const gsMatrix<T> &deriv2,
                                            const std::pair<short_t, short_t> & Dim,
                                            gsMatrix<T> &result)
    {
        const short_t domDim    = Dim.first;
        const int blockSize = ((domDim+1)*domDim)/2;
        const int numBlocks = deriv2.rows()/blockSize;
        result.resize(numBlocks,deriv2.cols());
        for (int b=0; b< numBlocks; ++b)
            result.row(b)=deriv2.middleRows(blockSize*b,domDim).colwise().sum();
    }

    template <typename input, typename output>
    static inline void deriv2ToHessianSingle(
            const input   secDers,
            output  hessian,
            short_t domDim)
    {
        switch ( domDim )
        {
        case 1:
            hessian(0,0)=secDers(0,0);
            break;
        case 2:
            hessian(0,0)=secDers(0,0);
            hessian(0,1)=secDers(2,0);
            hessian(1,0)=secDers(2,0);
            hessian(1,1)=secDers(1,0);
            break;
        case 3:
            hessian(0,0)=secDers(0,0);
            hessian(0,1)=secDers(3,0);
            hessian(0,2)=secDers(4,0);

            hessian(1,0)=secDers(3,0);
            hessian(1,1)=secDers(1,0);
            hessian(1,2)=secDers(5,0);

            hessian(2,0)=secDers(4,0);
            hessian(2,1)=secDers(5,0);
            hessian(2,2)=secDers(2,0);
            break;
        case 4:
            hessian(0,0)=secDers(0,0);
            hessian(0,1)=secDers(4,0);
            hessian(0,2)=secDers(5,0);
            hessian(0,3)=secDers(6,0);

            hessian(1,0)=secDers(4,0);
            hessian(1,1)=secDers(1,0);
            hessian(1,2)=secDers(7,0);
            hessian(1,3)=secDers(8,0);

            hessian(2,0)=secDers(5,0);
            hessian(2,1)=secDers(7,0);
            hessian(2,2)=secDers(2,0);
            hessian(2,3)=secDers(9,0);

            hessian(3,0)=secDers(6,0);
            hessian(3,1)=secDers(8,0);
            hessian(3,2)=secDers(9,0);
            hessian(3,3)=secDers(3,0);
            break;
        default:
            GISMO_ERROR("curl is not implemented in 5d or more");
        }
    }

    template <typename input>
    static inline void deriv2ToHessianSingle(
            const input         secDers,
            gsMatrix<T>  &hessian,
            short_t domDim)
    {
        switch ( domDim )
        {
        case 1:
            hessian(0,0)=secDers(0,0);
            break;
        case 2:
            hessian(0,0)=secDers(0,0);
            hessian(0,1)=secDers(2,0);
            hessian(1,0)=secDers(2,0);
            hessian(1,1)=secDers(1,0);
            break;
        case 3:
            hessian(0,0)=secDers(0,0);
            hessian(0,1)=secDers(3,0);
            hessian(0,2)=secDers(4,0);

            hessian(1,0)=secDers(3,0);
            hessian(1,1)=secDers(1,0);
            hessian(1,2)=secDers(5,0);

            hessian(2,0)=secDers(4,0);
            hessian(2,1)=secDers(5,0);
            hessian(2,2)=secDers(2,0);
            break;
        case 4:
            hessian(0,0)=secDers(0,0);
            hessian(0,1)=secDers(4,0);
            hessian(0,2)=secDers(5,0);
            hessian(0,3)=secDers(6,0);

            hessian(1,0)=secDers(4,0);
            hessian(1,1)=secDers(1,0);
            hessian(1,2)=secDers(7,0);
            hessian(1,3)=secDers(8,0);

            hessian(2,0)=secDers(5,0);
            hessian(2,1)=secDers(7,0);
            hessian(2,2)=secDers(2,0);
            hessian(2,3)=secDers(9,0);

            hessian(3,0)=secDers(6,0);
            hessian(3,1)=secDers(8,0);
            hessian(3,2)=secDers(9,0);
            hessian(3,3)=secDers(3,0);
            break;
        default:
            GISMO_ERROR("curl is not implemented in 5d or more");
        }
    }

    template <typename input, typename output>
    static inline void hessianToDeriv2Single(
            const input       hessian,
            output      deriv2,
            short_t domDim)
    {
        switch ( domDim )
        {
        case 1:
            deriv2(0,0)=hessian(0,0);
            break;
        case 2:
            deriv2(0,0)=hessian(0,0);
            deriv2(1,0)=hessian(1,1);
            deriv2(2,0)=(hessian(0,1)+hessian(1,0))/2;
            break;
        case 3:
            deriv2(0,0)=hessian(0,0);
            deriv2(1,0)=hessian(1,1);
            deriv2(2,0)=hessian(2,2);
            deriv2(3,0)=(hessian(0,1)+hessian(1,0))/2;
            deriv2(4,0)=(hessian(2,0)+hessian(2,0))/2;
            deriv2(5,0)=(hessian(1,2)+hessian(2,1))/2;
            break;
        case 4:
            deriv2(0,0)=hessian(0,0);
            deriv2(1,0)=hessian(1,1);
            deriv2(2,0)=hessian(2,2);
            deriv2(3,0)=hessian(3,3);

            deriv2(4,0)=(hessian(0,1)+hessian(1,0))/2;
            deriv2(5,0)=(hessian(2,0)+hessian(2,0))/2;
            deriv2(6,0)=(hessian(3,0)+hessian(3,0))/2;
            deriv2(7,0)=(hessian(1,2)+hessian(2,1))/2;
            deriv2(8,0)=(hessian(1,3)+hessian(3,1))/2;
            deriv2(8,0)=(hessian(2,3)+hessian(3,2))/2;
            break;
        default:
            GISMO_ERROR("curl is not implemented in 5d or more");
        }
    }


    template <typename input>
    static inline void hessianToDeriv2Single(
            const input  hessian,
            gsMatrix<T> &deriv2,
            short_t domDim)
    {
        switch ( domDim )
        {
        case 1:
            deriv2(0,0)=hessian(0,0);
            break;
        case 2:
            deriv2(0,0)=hessian(0,0);
            deriv2(1,0)=hessian(1,1);
            deriv2(2,0)=(hessian(0,1)+hessian(1,0))/2;
            break;
        case 3:
            deriv2(0,0)=hessian(0,0);
            deriv2(1,0)=hessian(1,1);
            deriv2(2,0)=hessian(2,2);
            deriv2(3,0)=(hessian(0,1)+hessian(1,0))/2;
            deriv2(4,0)=(hessian(2,0)+hessian(2,0))/2;
            deriv2(5,0)=(hessian(1,2)+hessian(2,1))/2;
            break;
        case 4:
            deriv2(0,0)=hessian(0,0);
            deriv2(1,0)=hessian(1,1);
            deriv2(2,0)=hessian(2,2);
            deriv2(3,0)=hessian(3,3);

            deriv2(4,0)=(hessian(0,1)+hessian(1,0))/2;
            deriv2(5,0)=(hessian(2,0)+hessian(2,0))/2;
            deriv2(6,0)=(hessian(3,0)+hessian(3,0))/2;
            deriv2(7,0)=(hessian(1,2)+hessian(2,1))/2;
            deriv2(8,0)=(hessian(1,3)+hessian(3,1))/2;
            deriv2(8,0)=(hessian(2,3)+hessian(3,2))/2;
            break;
        default:
            GISMO_ERROR("curl is not implemented in 5d or more");
        }
    }



    static inline  void deriv2ToHessian  (const gsMatrix<T> &deriv2,
                                          const std::pair<short_t, short_t> & Dim,
                                          gsMatrix<T> &result)
    {
        const short_t domDim=Dim.first;
        const int blockSize = ((domDim+1)*domDim)/2;
        const int numBlocks = deriv2.rows()/blockSize;
        result.resize(numBlocks*domDim,domDim*deriv2.cols());
        for (int p=0; p<deriv2.cols(); ++p)
            for (int b=0; b< numBlocks; ++b)
                deriv2ToHessianSingle(
                            deriv2.block(blockSize*b,p,blockSize,1),
                            result.block(domDim*b,domDim*p,domDim,domDim),
                            domDim
                            );
    }

    static inline  void  derivToMeasure(const gsMatrix<T> &deriv,
                                        const std::pair<short_t, short_t> & Dim,
                                        gsMatrix<T> &result)
    {
        //Note: suboptimal implementation for Dim.second==Dim.first
        const index_t blockSize = (Dim.first*Dim.second);

        GISMO_ASSERT( deriv.rows() % blockSize == 0,
                      "derivToMeasure: deriv dimensions do not agree");

        result.resize(deriv.rows()/blockSize, deriv.cols() );

        for (index_t b = 0; b < deriv.rows(); b+=blockSize) // for all blocks
            for (index_t p = 0; p != deriv.cols(); ++p) // for all points
            {
                // The transposed Jacobian matrix in the current column-block
                const gsAsConstMatrix<T> jacT(&deriv.coeffRef(b,p),
                                              Dim.first, Dim.second);
                result(b/blockSize,p) = math::sqrt( ( jacT * jacT.transpose() ).determinant() );
            }
    }

    template <short_t domainDim>
    static inline  void  dispatchPartialDerivToGradTransform(const gsMatrix<T> &deriv,
                                                             const std::pair<short_t, short_t> & Dim,
                                                             gsMatrix<T> &result)
    {
        switch (Dim.second)
        {
        case 4:  return derivToGradTransform<domainDim,4>(deriv,Dim,result);
        case 3:  return derivToGradTransform<domainDim,3>(deriv,Dim,result);
        case 2:  return derivToGradTransform<domainDim,2>(deriv,Dim,result);
        case 1:  return derivToGradTransform<domainDim,1>(deriv,Dim,result);
        default: return derivToGradTransform<domainDim,-1>(deriv,Dim,result);
        }
    }


    static inline  void  dispatchDerivToGradTransform(const gsMatrix<T> &deriv,
                                                      const std::pair<short_t, short_t> & Dim,
                                                      gsMatrix<T> &result)
    {
        switch (Dim.first)
        {
        case 4:  return dispatchPartialDerivToGradTransform<4>(deriv,Dim,result);
        case 3:  return dispatchPartialDerivToGradTransform<3>(deriv,Dim,result);
        case 2:  return dispatchPartialDerivToGradTransform<2>(deriv,Dim,result);
        case 1:  return dispatchPartialDerivToGradTransform<1>(deriv,Dim,result);
        default: return dispatchPartialDerivToGradTransform<-1>(deriv,Dim,result);
        }
    }

    template <short_t parDim, short_t tarDim>
    static inline  void  derivToGradTransform(const gsMatrix<T> &deriv,
                                              const std::pair<short_t, short_t> & Dim,
                                              gsMatrix<T> &result)
    {
        GISMO_ASSERT(parDim==-1||parDim==Dim.first, "wrong dimension");
        GISMO_ASSERT(tarDim==-1||tarDim==Dim.second, "wrong dimension");

        const int blockSize = (Dim.first*Dim.second);
        const short_t targetDim = Dim.second;
        const short_t domainDim = Dim.first;
        const int numBlocks = deriv.rows()/blockSize;
        result.resizeLike(deriv);
        for (index_t p=0; p<deriv.cols(); ++p)
            for (index_t b=0; b<numBlocks;++b)
            {
                Eigen::Map<const Eigen::Matrix<T,tarDim,parDim,Eigen::RowMajor> > jac(&deriv.coeffRef(b*blockSize,p),targetDim,domainDim);
                Eigen::Map<Eigen::Matrix<T,tarDim,parDim> > res(&result(b*blockSize,p),targetDim,domainDim);
                res = jac*(jac.transpose()*jac).inverse().eval();
            }
    }

    static inline  void derivToNormalSingle (Eigen::Map<const Eigen::Matrix<T,-1,-1,Eigen::RowMajor> > &jac, Eigen::Map<Eigen::Matrix<T,-1,-1> > &res)
    {
        GISMO_ASSERT(jac.rows()-jac.cols()==1, "not implemented");
        switch (jac.rows())
        {
        case 2:
            res(0,0)=-jac(0,1);
            res(1,0)=jac(0,0);
            break;
        case 3:
            // I hate Eigen!! jac.col(0).cross(jac.col(1)) is not compiling because they check that the size is fixed to 3 at compile time
            res(0,0)=jac(1,0)*jac(2,1)-jac(2,0)*jac(1,1);
            res(1,0)=jac(2,0)*jac(0,1)-jac(0,0)*jac(2,1);
            res(2,0)=jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);
            break;
        default:
            GISMO_ERROR("not implemented");
        }
    }


    static inline  void  derivToNormal(const gsMatrix<T> &deriv,
                                       const std::pair<short_t, short_t> & Dim,
                                       gsMatrix<T> &result)
    {
        const int inBlockSize  = (Dim.first*Dim.second);
        const short_t targetDim = Dim.second;
        const short_t domainDim = Dim.first;
        const int outBlockSize = targetDim*(targetDim-domainDim);

        GISMO_ASSERT(domainDim==1 || domainDim==2, "error, implemented only for curves and surfaces");
        GISMO_ASSERT(targetDim-domainDim==1, "error, implemented only in codim 1");

        const int numBlocks = deriv.rows()/inBlockSize;
        result.resize(outBlockSize*numBlocks,deriv.cols());

        for (index_t p=0; p<deriv.cols(); ++p)
            for (index_t b=0; b<numBlocks;++b)
            {
                Eigen::Map<const Eigen::Matrix<T,-1,-1,Eigen::RowMajor> > jac(&deriv.coeffRef(b*inBlockSize,p),targetDim,domainDim);
                Eigen::Map<Eigen::Matrix<T,-1,-1> > res(&result.coeffRef(b*inBlockSize,p),targetDim,domainDim);
                derivToNormalSingle(jac,res);
            }
    }

}; // convertValue struct

} // namespace gismo


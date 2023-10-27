/** @file gsDomainDecomposition.h

    @brief Provides multi-patch preconditioners, particularly smoothers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsPreconditioner.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsCore/gsMultiPatch.h>

namespace gismo
{

/// gsBoundary
struct boxComponent2 {
    index_t m_index;
    index_t m_total_dim;

    boxComponent2( index_t b, index_t total_dim ) : m_index(b), m_total_dim(total_dim) {}
    boxComponent2( const boxSide& p, index_t total_dim ) : m_total_dim(total_dim)
    {
        const index_t d = (p.m_index-1)/2;
        const index_t o = (p.m_index-1)%2;
        //m_index = (o+1)*3^d
        m_index = o+1;
        for (index_t c = 0; c<d; ++c)
            m_index *= 3;
    }
    boxComponent2( boxCorner p, index_t total_dim ) : m_total_dim(total_dim)
    {
        p.m_index -= 1;
        m_index = 0;
        index_t c = 1;
        for (index_t d=0; d<total_dim; ++d)
        {
            m_index += (1+p.m_index%2) * c;
            c *= 3;
        }
    }
    index_t dim() const
    {
         index_t result = 0, tmp = m_index;
         for (index_t i=0; i<m_total_dim; ++i)
         {
             if (tmp%3 == 0)
                 ++result;
             tmp /= 3;
         }
         return result;
    }
    std::vector<boxCorner> containedCorners() const
    {
        const index_t result_sz = 1u<<dim();

        std::vector<boxCorner> result;
        result.reserve(result_sz);

        for (index_t i=0; i<result_sz; ++i)
        {
             index_t idx = 0;
             index_t ii=i, jj=m_index, c=1;
             for (index_t d=0; d<m_total_dim; ++d)
             {
                 if (jj%3==0)
                 {
                     idx += c*(ii%2);
                     ii /= 2; jj /= 3;
                 }
                 else
                 {
                     idx += c*(jj%3-1);
                     jj /= 3;
                 }
                 c *= 2;
             }
             result.push_back(boxCorner(1+idx));
        }
        return result;
    }

    boxSide asSide()
    {
        GISMO_ENSURE( dim() == m_total_dim-1, "This is not a side." );
        index_t d = 0, idx = m_index;
        while (idx>0)
        {
            if (idx%3)
                return boxSide( idx%3 + 2*d );
            idx /= 3;
            ++d;
        }
        throw std::runtime_error("idx was 0");
    }

    boxCorner asCorner()
    {
        GISMO_ENSURE( dim() == 0, "This is not a corner." );
        index_t factor = 1, idx = m_index, result = 0;
        while (idx>0)
        {
            result += ( idx%3 - 1 ) * factor;
            idx /= 3;
            factor *= 3;
        }
        return boxCorner(1+result);
    }

    enum location {
        interior = 0,
        begin = 1,
        end = 2
    };

    location parameter(index_t direction)
    {
        GISMO_ASSERT( direction >= 0 && direction < m_total_dim, "Out of bounds" );
        index_t result = m_index;
        for (index_t i=0; i<direction; ++i)
            result /= 3;
        return location(result%3);
    }

    void setParameter(index_t direction, location par)
    {
        const index_t diff = par - parameter(direction);
        if (diff)
        {
            index_t factor = 1;
            for (index_t i=0; i<direction; ++i)
                factor *= 3;
            m_index += diff * factor;
        }
    }

};

struct patchComponent2 : boxComponent2 {
    index_t patch;
    patchComponent2( index_t p, index_t b, index_t total_dim ) : boxComponent2(b,total_dim), patch(p) {}
    patchComponent2( index_t p, boxComponent2 b ) : boxComponent2(b), patch(p) {}
    patchComponent2( const patchSide& p, index_t total_dim ) : boxComponent2(p,total_dim), patch(p) {}
    patchComponent2( patchCorner p, index_t total_dim ) : boxComponent2(p, total_dim), patch(p) {}
    std::vector<patchCorner> containedCorners() const
    {
        std::vector<boxCorner> tmp = boxComponent2::containedCorners();
        std::vector<patchCorner> result;
        const index_t sz = tmp.size();
        result.reserve(sz);
        for (index_t i=0; i<sz; ++i)
            result.push_back(patchCorner(patch,tmp[i]));
        return result;
    }
    patchSide asSide()
    {
        return patchSide( patch, boxComponent2::asSide() );
    }
    patchCorner asCorner()
    {
        return patchCorner( patch, boxComponent2::asCorner() );
    }
};

///


/// gsBasis

template<class T>
struct gsBasisHelper {

    typedef typename gsBasis<T>::uPtr   BasisUPtr;
    gsMatrix<index_t> boundary(const gsBasis<T>& basis, boxComponent2 b, bool strict = false);
    BasisUPtr boundaryBasis(const gsBasis<T>& basis, boxComponent2 b);

    static gsSparseMatrix<T,RowMajor> setupTransferMatrix( const gsMatrix<index_t>& coefs, index_t totalSize );
    static gsSparseMatrix<T,RowMajor> setupTransferMatrix( const std::vector<index_t>& coefs, index_t totalSize );


};

/// gsMultiBasis

template<class T>
struct gsMultiBasisHelper {

    struct component {
        std::vector<patchComponent2> components;
        std::vector<patchCorner> corners;
        void swap( component& other ) { components.swap(other.components); corners.swap(other.corners); }
    };

    static std::vector<component>
    getComponents( const gsMultiBasis<T>& mb, bool combineCorners = false );

    static gsSparseMatrix<T, RowMajor>
    getTransfer(
        const std::vector<patchComponent2>& pc,
        const gsMultiBasis<T>& mb,
        const gsBoundaryConditions<T>& bc,
        const gsOptionList& opt,
        bool strict = false
    );

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsDomainDecomposition.hpp)
#endif

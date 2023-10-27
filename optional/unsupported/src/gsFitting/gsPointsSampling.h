/** @file gsPointsSampling.h

    @brief Contains the class gsPointsSampling that permits to samples
    points of a given domain using a given basis.

    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include "gsFittingIterator.h"

namespace gismo
{


/**
   @brief class gsPointsSampling: sample the domain of a gsBasis
   using the Gauss quadrature points. The advantage here is that
   the sampling is done according to the degrees of freedom of the basis.
   On the other side, before using this class to generate a sampling,
   the basis must be refined for patches with very few degrees of freedom.
   The basis can be defined on a single patch (gsBasis) or on a
   multipatch (gsMultiBasis).
   \ingroup
**/

template<class GenBasis, class T>
class gsPointsSampling
{

public:
    /// Constructor: sets the basis and the local size.
    /// local_size: the number of points chosen
    /// in each elementary cell for the sampling
    gsPointsSampling(GenBasis &basis, int local_size) :
    m_basis(basis), m_numNodes(basis.domainDim()),
    m_quRule((init_numNodes(local_size), m_numNodes))
    {   init_gsPointsSampling();  }

    gsPointsSampling(GenBasis &basis, int* local_size = NULL) :
    m_basis(basis), m_numNodes(basis.domainDim()),
    m_quRule((init_numNodes(local_size), m_numNodes))
    {   init_gsPointsSampling();  }

    /// Destructor
    ~gsPointsSampling(){  }

    /// Performs the sampling using a domain iterator on the domain
    /// of the basis m_basis. The sampling is returned in pts.
    void points_sampling(std::vector< gsMatrix<T> >& pts);

    /// Returns the number of points to be sampled (or already sampled)
    /// in each patches
    std::vector<index_t>& sizes(){ return m_size_sample;  }

    /// Returns the number of points to be sampled (or already sampled)
    /// in some patch
    index_t size(int patch){ return m_size_sample[patch];  }

    /// Returns the total number of points to be sampled (or already sampled)
    index_t total_size();

private:
    /// Initialize the class, notably the number of points to be sampled.
    void init_gsPointsSampling();

    /// Sets the number of points to be sampled in each elementary cell.
    void init_numNodes(int local_size);
    void init_numNodes(int* local_size);

    /// The basis used for the sampling
    GenBasis &m_basis;

    /// The number of points to be sampled
    /// (or already sampled) in each patches
    std::vector<index_t> m_size_sample;

    /// Sets the number of points to be sampled in each elementary cell
    /// (for each dimension)
    gsVector<int> m_numNodes;

    /// The Gaussian rule used for the sampling
    gsGaussRule<T> m_quRule;
};


}// namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPointsSampling.hpp)
#endif

//////////////////////////////////////////////////
//////////////////////////////////////////////////

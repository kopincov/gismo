/** @file gsPointsSampling.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsFitting/gsPointsSampling.h>
#include <gsFitting/gsFittingUtilsGen.h>

namespace gismo
{

template<class GenBasis, class T>
void gsPointsSampling<GenBasis, T>::
points_sampling(std::vector< gsMatrix<T> >& pts)
{
    gsGenDomainIterator<T> domIt(m_basis);

    unsigned dim = domIt.paramDim();
    unsigned nbpts_loc = m_quRule.numNodes();

    gsMatrix<T> quNodes(dim, nbpts_loc);
    gsVector<T> weights_unused(nbpts_loc);
    int patch = 0, patch_prec = 0;
    int ind = 0;

    for (; domIt.good(); domIt.next() )
    {
        m_quRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                        quNodes, weights_unused);
        patch = domIt.patch();
        if(patch != patch_prec)
            ind = 0;
        for(unsigned i = 0;i < nbpts_loc;i++)
        {
            pts[patch].col(ind) = quNodes.col(i);
            ind++;
        }

        patch_prec = patch;
    }
}

template<class GenBasis, class T>
index_t gsPointsSampling<GenBasis, T>::
total_size()
{
    unsigned res = 0;
    unsigned s = m_size_sample.size();;
    for(unsigned i = 0;i < s;i++)
        res += m_size_sample[i];
    return res;
}


template<class GenBasis, class T>
void gsPointsSampling<GenBasis, T>::
init_gsPointsSampling()
{
    gsGenDomainIterator<T> tmp(m_basis);
    unsigned s = tmp.nPatches();
    unsigned nbpts_loc = m_quRule.numNodes();
    for(unsigned i = 0;i < s;i++)
    {
        unsigned nb_elements = tmp.numElement(i);
        m_size_sample.push_back(nbpts_loc * nb_elements);
    }
}

template<class GenBasis, class T>
void gsPointsSampling<GenBasis, T>::
init_numNodes(int local_size)
{
    unsigned dim = m_basis.domainDim();
    for (unsigned i = 0;i != dim;i++)
    {
        if(local_size == -1)
            m_numNodes[i] = firstBasis(m_basis).degree(i);
        else
            m_numNodes[i] = local_size;
    }
}

template<class GenBasis, class T>
void gsPointsSampling<GenBasis, T>::
init_numNodes(int* local_size)
{
    unsigned dim = m_basis.domainDim();
    for (unsigned i = 0;i != dim;i++)
    {
        if(local_size == NULL)
            m_numNodes[i] = firstBasis(m_basis).degree(i);
        else
            m_numNodes[i] = local_size[i];
    }
}

} // namespace gismo

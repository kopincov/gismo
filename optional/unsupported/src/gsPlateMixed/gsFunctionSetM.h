/** @file gsFunctionSet.h

    @brief Function to compute the bending moments M

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/

#pragma once

#include <gsCore/gsFunctionSet.h>
#include <gsPlateMixed/gsFunctionM.h>

namespace gismo {

template <typename T>
class gsFunctionSetM : public gsFunctionSet<T>
{
public:
    /// Shared pointer for gsFunctionSetM
    typedef memory::shared_ptr< gsFunctionSetM > Ptr;

    /// Unique pointer for gsFunctionSetM
    typedef memory::unique_ptr< gsFunctionSetM > uPtr;

public:

    gsFunctionSetM(gsMultiPatch<T>& mp, gsMultiPatch<T>& p, gsMultiPatch<T>& phi)
    {
        M.resize(mp.size());
        for(index_t i= 0; i!=mp.size(); ++i)
        {
            M[i] = memory::make_shared(new gsFunctionM<T>(mp.piece(i), p.piece(i), phi.piece(i)));
        }

    }

    static uPtr make(gsMultiPatch<T>& mp, gsMultiPatch<T>& u_co) {return memory::make_unique(new gsFunctionSetM(mp,u_co)); }

    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsFunctionSetM, clone)

    /// @brief Returns the piece(s) of the function(s) at subdomain \a k
    virtual const gsFunctionM<T> & piece(const index_t k) const
    {
        return *M[k];
    }


    /**
       @brief Dimension of the (source) domain.
       @return For \f$f:\mathbb{R}^n\rightarrow\mathbb{R}^m\f$, returns \f$n\f$.
     */
    virtual short_t domainDim () const
    {
        return 2;
    }

private:

    std::vector<typename gsFunctionM<T>::Ptr > M;


};



} // namespace gismo



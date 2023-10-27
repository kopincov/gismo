/** @file gsFunctionSet.h

    @brief Function to compute the physical displacement

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): K. Rafetseder
*/

#pragma once

#include <gsCore/gsFunctionSet.h>
#include <gsPlateMixed/gsFunctionPDisplace.h>

namespace gismo {

template <typename T>
class gsFunctionSetPDisplace : public gsFunctionSet<T>
{
public:
    /// Shared pointer for gsFunctionSetPDisplace
    typedef memory::shared_ptr< gsFunctionSetPDisplace > Ptr;

    /// Unique pointer for gsFunctionSetPDisplace
    typedef memory::unique_ptr< gsFunctionSetPDisplace > uPtr;

public:

    gsFunctionSetPDisplace(gsMultiPatch<T>& mp, gsMultiPatch<T>& u_co)
    {
        uPhysical.resize(mp.nPatches());
        for(size_t i= 0; i!=mp.nPatches(); ++i)
        {
            uPhysical[i] = memory::make_shared(new gsFunctionPDisplace<T>(mp.piece(i), u_co.piece(i)));
        }

    }

    static uPtr make(gsMultiPatch<T>& mp, gsMultiPatch<T>& u_co) {return memory::make_unique(new gsFunctionSetPDisplace(mp,u_co)); }

    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsFunctionSetPDisplace, clone)

    /// @brief Returns the piece(s) of the function(s) at subdomain \a k
    virtual const gsFunctionPDisplace<T> & piece(const index_t k) const
    {
        return *uPhysical[k];
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

    std::vector<typename gsFunctionPDisplace<T>::Ptr > uPhysical;


};



} // namespace gismo



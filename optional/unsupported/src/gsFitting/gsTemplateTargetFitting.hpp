/** @file gsTemplateTargetFitting.hpp

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#include <gsFitting/gsTemplateTargetFitting.h>

namespace gismo
{

template<class T> gsFunctionSet<T>*
gsTemplateTargetFitting<T>::
computeMappingTrimmingDyn(gsMultiBasis<T>* basis)
{
    index_t dim = dim_dom();
    if(dim == 2)
        return computeMappingTrimming<2>(basis);
    else
    {
        GISMO_ENSURE(dim == 3, "We only consider dimensions 2 and 3");
        return computeMappingTrimming<3>(basis);
    }
}

template<class T> template<unsigned d>
gsFunctionSet<T>* gsTemplateTargetFitting<T>::
computeMappingTrimming(gsMultiBasis<T>* basis)
{
    index_t _dim_im = dim_im();
    if(basis == NULL)
    {
        gsBasis<T>* p_basis;
        gsTensorBSplineBasis<d, T> basis2
            = m_coupling.template
            getBasisTrimming<d>(*m_param_fitting);
        gsTHBSplineBasis<d> thb(basis2);

        if(m_param_fitting->use_refinement)
            p_basis = &thb;
        else
            p_basis = &basis2;

        return fitData<d, gsBasis<T>, T>
            (_dim_im, *p_basis, *m_param_fitting, m_coupling);
    }
    else
    {
        GISMO_ENSURE(basis->nBases() == (unsigned) m_template.nPatches(),
                     "The number of patches does not correspond. The bc should be split maybe? ");
        return fitData<d, gsMultiBasis<T>, T>
            (_dim_im, *basis, *m_param_fitting, m_coupling);
    }

}

} /// namespace gismo

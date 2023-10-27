/** @file gsFittingIdConstr.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsParameterLines.h>
#include <gsFitting/gsFittingUtilsCheck.hpp>
#include <gsFitting/gsFittingIdConstr.h>


namespace gismo
{

template<short_t d, class GenBasis, class T>
_gsFittingIdConstr<d, GenBasis, T>::
_gsFittingIdConstr(bool export_iter, std::string& output,
                          int degree)
: m_output(output)
{
    m_export = export_iter;
    m_degree = degree;
    m_template = NULL;
}


template<short_t d, class GenBasis, class T>
gsFittingIdConstr<d, GenBasis, T>::
gsFittingIdConstr(GenBasis &model, bool export_iter,
                  std::string& output,
                  int degree) :
_gsFittingIdConstr<d, GenBasis, T>::
_gsFittingIdConstr(export_iter, output, degree){ }


template<short_t d, class T>
gsFittingIdConstr<d, gsBasis<T>, T>::
gsFittingIdConstr(gsBasis<T> &model, bool export_iter,
                  std::string& output,
                  int degree) :
_gsFittingIdConstr<d, gsBasis<T>, T>::
_gsFittingIdConstr(export_iter, output, degree)
{
    init_basis(model);
}


template<short_t d, class T>
gsFittingIdConstr<d, gsMultiBasis<T>, T>::
gsFittingIdConstr(gsMultiBasis<T> &model, bool export_iter,
                  std::string& output, int degree) :
_gsFittingIdConstr<d, gsMultiBasis<T>, T>::
_gsFittingIdConstr(export_iter, output, degree)
{
    init_basis(model);
}

/// We have to consider separatly the single patch case
/// and the multipatch case for the construction of the basis
template<short_t d, class T> void
gsFittingIdConstr<d, gsBasis<T>, T>::
init_basis(gsBasis<T> &model)
{
    /// the knot vectors of the template
    std::vector< gsKnotVector<T> > m_kv;

    gsMatrix<T> support = model.support();
    unsigned s = support.rows();
    for(unsigned i = 0;i < s;i++)
    {
        m_kv.push_back(
            gsKnotVector<T>(support(i,0), support(i,1),
                            0, m_degree + 1) );
    }
    m_basis = gsTensorBSplineBasis<d>(m_kv);
}


template<short_t d, class T> void
gsFittingIdConstr<d, gsMultiBasis<T>, T>::
init_basis(gsMultiBasis<T> &model)
{
    const gsBoxTopology& topology ( model.topology() );

    /// the knot vectors of the template
    unsigned nb_patches = model.nBases();
    std::vector< gsKnotVector<T> > m_kv;
    std::vector<gsBasis<T>*> vect_point;

    for(unsigned i = 0;i < nb_patches;i++)
    {
        gsMatrix<T> support = model[i].support();
        unsigned s = support.rows();
        for(unsigned i = 0;i < s;i++)
        {
            m_kv.push_back(
                gsKnotVector<T>(support(i,0), support(i,1),
                                0, m_degree + 1) );
        }
        vect_point.push_back(new gsTensorBSplineBasis<d, T>(m_kv));
        m_kv.clear();
    }
    m_basis = gsMultiBasis<T>(vect_point, topology);
}


/// We suppose for the moment that the template domain is
/// given and we take the identity in the other case
template<short_t d, class GenBasis, class T>
void _gsFittingIdConstr<d, GenBasis, T>::
init_identity_mapping(unsigned dim_im,
                      gsFittingEnergy<GenBasis, T>& ener)
{
    gsFittingIdConstr<d, GenBasis>* _this
        = static_cast<gsFittingIdConstr
                      <d, GenBasis, T>* >(this);
    GenBasis &basis ( _this->basis() );
    gsFittingParam<T> param;
    bool remove_result = false;
    param.deform_min = false;

    gsFittingBase<d, GenBasis, T>
        fitt_id(dim_im, basis, remove_result, param);
    ener.copy_non_smoothing(fitt_id.energy());

    m_template = compute(fitt_id);

}


template<short_t d, class GenBasis, class T>
gsFunctionSet<T>* _gsFittingIdConstr<d, GenBasis, T>::
compute(gsFittingBase<d, GenBasis, T>& fitt)
{
    GenBasis& basis = fitt.basis();
    int dim = fitt.dim_im();

    bool finished = false;
    int ind = 0;
    gsFittingIntegrandLin<d, 2, T>* integr = NULL;
    T coeff_smoo = 1.;
    gsFunctionSet<T>* result = NULL;
    int nb_max = 8;

    ///// DBG //////
    m_export = false;
    ////////////////
    while(! finished)
    {
        fitt.compute();
        /// we check that the determinant is positive at each corner.
        /// If not, we recompute the template
        result = fitt.result();
        GISMO_ASSERT(result != NULL, "Error during the computation of the initialization of the iterative algorithm");
        if(positiveDet<T>(basis, *result))
            finished = true;
        if(! finished)
        {
            std::string parameter_lines_filename
                = m_output + "_template_iter_"
                + util::to_string(ind);
            if(m_export || ind > nb_max)
            {
                writeParameterLines(getGeometry(basis, *result),
                                    parameter_lines_filename,
                                    50, 3);
            }
            if(ind > nb_max)
            {
                GISMO_ERROR("Error during the construction of the initialization. Maybe, the template and the target do not match (different orientation,...). You can check on the file: "
                            + parameter_lines_filename);
                return NULL;
            }
            fitt.reset();

            if(ind == 0)
            {
                std::vector< gsFittingQuadrature<d, GenBasis,
                                                 2, T> >&
                    vect_quadr(fitt.quadrEner().template
                               getQuadratures<d>());
                GISMO_ASSERT(vect_quadr.size() == 0,
                             "This vector should be empty");
                vect_quadr.push_back(gsFittingQuadrature<d, GenBasis, 2, T>(basis, dim));
                integr = new gsFittingIntegrandLin<d, 2, T>
                    (coeff_smoo, 0., 1., true);
                vect_quadr[0].addIntegrand(integr);
            }
            else
            {
                coeff_smoo *= 10.;
                integr->setCoeffGlobal(coeff_smoo);
            }
            ind++;
        }
    }
    return fitt.result();
}

} // namespace gismo

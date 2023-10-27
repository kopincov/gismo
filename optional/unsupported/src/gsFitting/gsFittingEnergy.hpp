/** @file gsFittingEnergy.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/
#pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <gsFitting/gsFittingEnergy.h>
#include <gsFitting/gsParameterLines.h>

#include <gsFitting/gsFittingQuadrature.hpp>
#include <gsFitting/gsLeastSquares.hpp>
#include <gsFitting/gsLeastSquaresLMult.hpp>
#include <gsFitting/gsFittingUtilsGen.h>
#include <gsFitting/gsFittingUtilsGen.hpp>

#include <gsFitting/gsFittingSystem.hpp>

namespace gismo
{

template<class GenBasis, class T>
void gsFittingQuadrEnergy<GenBasis, T>::
set_properties_integrands()
{
    m_split_dim = true;
    m_isLinear = true;

    index_t s;
    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(index_t i = 0;i < s;i++)
        {
            gsFittingQuadratureGen<GenBasis, 2, T>&
                quadr = getQuadrature(d, i);
            m_split_dim = m_split_dim
                && quadr.canSplitDimension();
            m_isLinear = m_isLinear && quadr.isLinear();
        }
    }
}

template<class GenBasis, class T>
bool gsFittingQuadrEnergy<GenBasis, T>::has_smoothing()
{
    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            if(getQuadrature(d, i).has_smoothing())
                return true;
    }
    return false;
}


template<class GenBasis, class T>
void gsFittingQuadrEnergy<GenBasis, T>::
assemble(gsFittingSystem<GenBasis, T>& system)
{
    unsigned s;
    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            getQuadrature(d, i).
                assembleSystem(system, m_isLinear);
    }
}


template<class GenBasis, class T>
void gsFittingQuadrEnergy<GenBasis, T>::
set_current_map(gsFunctionSet<T>* current)
{
    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            getQuadrature(d, i).set_current_map(current);
    }
}

template<class GenBasis, class T>
void gsFittingQuadrEnergy<GenBasis, T>::
decrease_coeff_smoothing(T ratio)
{
    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            getQuadrature(d, i).
                decrease_coeff_smoothing(ratio);
    }
}

template<class GenBasis, class T>
void gsFittingQuadrEnergy<GenBasis, T>::
decrease_linear_coeff_smoothing(T ratio)
{
    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            getQuadrature(d, i)
                .decrease_linear_coeff_smoothing(ratio);
    }
}


template<class GenBasis, class T> template<unsigned d>
void gsFittingQuadrEnergy<GenBasis, T>::
addGlobalIntegrand(gsFittingIntegrand<d,T>* integr,
                   index_t dim_im, GenBasis& basis)
{
    typename ContainerQuadr<d>::type&
        quadr = getQuadratures<d>();
    index_t s = quadr.size();
    if(s == 0)
    {
        quadr.push_back(gsFittingQuadrature<d, GenBasis, 2, T>(basis, dim_im));
        quadr[0].addIntegrand(integr);
    }
    else if(s == 1)
        quadr[0].addIntegrand(integr);

    else
        gsWarn << "There exists more than one global integrand. The integrand has not been added!!"
               << std::endl;
}


template<class GenBasis, class T> template<unsigned d>
void gsFittingQuadrEnergy<GenBasis, T>::
add_regularization_term(T coeff, index_t dim, GenBasis& basis)
{
    addGlobalIntegrand<d>
        ( new gsFittingIntegrandL2Dist<d, gsGeometry<T>, T>
          (dim, NULL, false, coeff), dim, basis);
    addGlobalIntegrand<d>
        ( new gsFittingIntegrandLin<d, 2, T>
          (coeff, 1., 0., false), dim, basis );
}

/// <------------- Energy computation --------------->


template<class GenBasis, class T>
void gsFittingEnergy<GenBasis, T>::
print_error(bool compute)
{
    if(compute)
        computeErrors();

    if(m_print_messages)
    {
        gsInfo << "-------- COMPUTE ERROR ------"
               << std::endl;
        gsInfo << "ERROR : max : " << maxError() << "   min : "
               << minError() << "    total : " << totalError()
               << "   (! norm l2 minimized and norm l1 printed)"
               << std::endl;

        if(has_smoothing())
        {
            gsInfo << "ENERGY : max : "
                   << m_quadrEner.maxEnergySmoothing()
                   << "   min : "
                   << m_quadrEner.minEnergySmoothing()
                   << "    total : "
                   << m_quadrEner.totalEnergySmoothing()
                   << std::endl;
        }
    }
}

template<class GenBasis, class T>
T gsFittingQuadrEnergy<GenBasis, T>::maxError()
{
    T res = 0.;
    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            res = std::max(res, getQuadrature(d, i).maxError());
    }
    return res;
}

template<class GenBasis, class T>
T gsFittingQuadrEnergy<GenBasis, T>::minError()
{
    T res = std::numeric_limits<T>::max();
    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            res = std::min(res, getQuadrature(d, i).minError());
    }
    return res;
}

template<class GenBasis, class T>
T gsFittingQuadrEnergy<GenBasis, T>::totalError()
{
    T res = 0;

    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            res += getQuadrature(d, i).totalError();
    }
    return res;
}


template<class GenBasis, class T>
T gsFittingQuadrEnergy<GenBasis, T>::
maxEnergySmoothing()
{
    T res = 0;
    T val;

    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
        {
            val = getQuadrature(d, i).maxEnergySmoothing();
            if(val < 0)
                return val;
            res += val;
        }
    }
    return res;
}

template<class GenBasis, class T>
T gsFittingQuadrEnergy<GenBasis, T>::
minEnergySmoothing()
{
    T res = 0.;
    T val;

    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
        {
            val = getQuadrature(d, i).minEnergySmoothing();
            if(val < 0.)
                return val;
            res += val;
        }
    }
    return res;
}

template<class GenBasis, class T>
T gsFittingQuadrEnergy<GenBasis, T>::totalEnergySmoothing()
{
    T res = 0;
    T val;

    unsigned s;

    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
        {
            val = getQuadrature(d, i).totalEnergySmoothing();
            if(val < 0)
                return val;
            res += val;
        }
    }
    return res;
}

template<class GenBasis, class T>
void gsFittingQuadrEnergy<GenBasis, T>::
computeErrors(index_t type)
{
    unsigned s;
    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            getQuadrature(d, i).computeEnergies(type);
    }
}

template<class GenBasis, class T>
void gsFittingEnergy<GenBasis, T>::
actualize_basis()
{
    m_least_squares.actualize_degrees_freedomLS();
}

template<class GenBasis, class T>
void gsFittingQuadrEnergy<GenBasis, T>::
actualize_identity_mapping(gsFunctionSet<T>* ident)
{
    unsigned s;
    for(unsigned d = 1;d <= 3;d++)
    {
        s = sizeQuadrature(d);
        for(unsigned i = 0;i < s;i++)
            getQuadrature(d, i).set_identity_mapping(ident);
    }
}

template<class GenBasis, class T>
void gsFittingQuadrEnergy<GenBasis, T>::
copy_non_smoothing(gsFittingQuadrEnergy<GenBasis, T>& new_ener)
{
    copy_non_smoothing<1>( getQuadratures<1>(),
                           new_ener.getQuadratures<1>() );
    copy_non_smoothing<2>( getQuadratures<2>(),
                           new_ener.getQuadratures<2>() );
    copy_non_smoothing<3>( getQuadratures<3>(),
                           new_ener.getQuadratures<3>() );
}



template<class GenBasis, class T> template<unsigned d>
void gsFittingQuadrEnergy<GenBasis, T>::
copy_non_smoothing(typename ContainerQuadr<d>::type& origin,
                   typename ContainerQuadr<d>::type& copy)
{
    unsigned s = origin.size();
    GISMO_ASSERT(copy.size() == 0,
                 "The container copy should be empty!!");
    index_t ind = 0;
    for(unsigned i = 0;i < s;i++)
    {
        if(origin[i].has_fitting())
        {
            copy.push_back(gsFittingQuadrature<d, GenBasis, 2, T>(origin[i]));
            copy[ind].removeAll();
            gsFittingQuadrature<d, GenBasis, 2, T>::
                copy_non_smoothing_integrand(origin[i], copy[ind]);
            ind++;
        }
    }
}

} // namespace gismo

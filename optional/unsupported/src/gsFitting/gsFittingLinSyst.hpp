/** @file gsFittingLinSyst.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/
#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsParameterLines.h>

#include <gsFitting/gsFittingLinSyst.h>
#include <gsFitting/gsFittingQuadrature.hpp>
#include <gsFitting/gsFittingIdConstr.hpp>
#include <gsFitting/gsLeastSquares.hpp>
#include <gsFitting/gsLeastSquaresLMult.hpp>
#include <gsFitting/gsFittingUtilsGen.hpp>


namespace gismo
{

template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::construct_template_mapping()
{
    gsFittingQuadrature<1, GenBasis, 2, T>* quadr_1d = NULL;
    gsFittingQuadrature<2, GenBasis, 2, T>* quadr_2d = NULL;
    gsFittingQuadrature<3, GenBasis, 2, T>* quadr_3d = NULL;

    GenBasis& basis_template ( m_template_constr.basis() );
    gsLeastSquares<GenBasis, T>* p_least_squares = NULL;
    if(m_least_squares != NULL)
    {
        p_least_squares = new gsLeastSquares<GenBasis, T>
            (basis_template, *m_least_squares, m_dim, m_output, false);
    }

    /// Copy the non smoothing integrands
    if(m_quadr_1d != NULL)
    {
        quadr_1d = static_cast
            <gsFittingQuadrature<1, GenBasis, 2, T>*>
            (m_quadr_1d->copy(&basis_template));
        bool active = gsFittingQuadrature<1, GenBasis, 2, T>::
            copy_non_smoothing_integrand(*m_quadr_1d, *quadr_1d);
        if(! active)
        {
            delete quadr_1d;
            quadr_1d = NULL;
        }
    }
    if(m_quadr_2d != NULL)
    {
        quadr_2d = static_cast
            <gsFittingQuadrature<2, GenBasis, 2, T>*>
            (m_quadr_2d->copy(&basis_template));
        bool active = gsFittingQuadrature<2, GenBasis, 2, T>::
            copy_non_smoothing_integrand(*m_quadr_2d, *quadr_2d);
        if(! active)
        {
            delete quadr_2d;
            quadr_2d = NULL;
        }
    }
    if(m_quadr_3d != NULL)
    {
        quadr_3d = static_cast
            <gsFittingQuadrature<3, GenBasis, 2, T>*>
            (m_quadr_3d->copy(&basis_template));
        bool active = gsFittingQuadrature<3, GenBasis, 2, T>::
            copy_non_smoothing_integrand(*m_quadr_3d, *quadr_3d);
        if(! active)
        {
            delete quadr_3d;
            quadr_3d = NULL;
        }
    }
    m_template_constr.init_identity_mapping(m_dim, p_least_squares,
                                            quadr_1d, quadr_2d,
                                            quadr_3d);
    if(quadr_1d != NULL)
        delete quadr_1d;
    if(quadr_2d != NULL)
        delete quadr_2d;
    if(quadr_3d != NULL)
        delete quadr_3d;
    if(p_least_squares != NULL)
        delete p_least_squares;
    m_template = m_template_constr.get_template();

}

template<short_t d, class GenBasis, class T> void
_gsFittingLinSyst<d, GenBasis, T>::
init_gsFittingLinSyst(GenBasis &basis, bool remove_result,
                  unsigned dim_im, gsFittingParam<T>& param)
 {
     m_remove_result = remove_result;
     m_result_linear = NULL;
     m_current = NULL;
     m_basis = &basis;
     m_deform_minim = param.deform_min;
     m_quadr_1d = NULL;
     m_quadr_2d = NULL;
     m_quadr_3d = NULL;
     m_least_squares = NULL;
     m_dim = dim_im;
     m_template = param.p_template;
     m_print_messages = param.print_messages;
     m_output = param.output;
     m_is_displacement = false;
     m_energy_max = 0.;
 }

template<short_t d, class GenBasis, class T>
_gsFittingLinSyst<d, GenBasis, T>::
_gsFittingLinSyst(unsigned dim_im, GenBasis &basis,
              bool remove_result,
              gsFittingParam<T>& param,
              gsLeastSquares<GenBasis, T>* least_squares)
: m_local_global(basis, dim_im), m_template_constr(basis, false,
                                                   param.output, 1),
  m_A_mat()
{
    init_gsFittingLinSyst(basis, remove_result, dim_im, param);
    m_least_squares = least_squares;

    if(m_template == NULL)
        m_template_given = false;
}

template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::set_properties_integrands()
{
    m_split_dim = true;
    m_isLinear = true;
    if(m_quadr_1d != NULL)
    {
        m_split_dim = m_split_dim && m_quadr_1d->canSplitDimension();
        m_isLinear = m_isLinear && m_quadr_1d->isLinear();
    }
    if(m_quadr_2d != NULL)
    {
        m_split_dim = m_split_dim && m_quadr_2d->canSplitDimension();
        m_isLinear = m_isLinear && m_quadr_2d->isLinear();
    }
    if(m_quadr_3d != NULL)
    {
        m_split_dim = m_split_dim && m_quadr_3d->canSplitDimension();
        m_isLinear = m_isLinear && m_quadr_3d->isLinear();
    }
    if(!m_isLinear)
    {
        if(m_template == NULL)
            m_template_given = false;
    }
}


template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::compute()
{
    this->set_properties_integrands();
    GISMO_ASSERT(m_result_linear == NULL, "Result already computed");
    if(m_deform_minim && m_template == NULL)
    {
        construct_template_mapping();
    }
    computeLinear();
}

template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::computeLinear()
{
    gsFittingLinSyst<d, GenBasis, T>* _this = static_cast
        <gsFittingLinSyst<d, GenBasis, T>*>(this);
    bool can_split = m_split_dim;

    int num_basis = sizeDof();
    if(! can_split)
        num_basis *= m_dim;

    //left side matrix
    gsSparseEntries<T> entries_mat;

    //right side vector (more dimensional!)
    gsMatrix<T> m_B;
    if(can_split)
        m_B = gsMatrix<T>(num_basis, m_dim);
    else
        m_B = gsMatrix<T>(num_basis, 1);
    // enusure that all entries are zero in the beginning
    m_B.setZero();
    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    if(m_least_squares != NULL)
        m_least_squares->assembleSystem(entries_mat, m_B, m_split_dim);

    if(m_print_messages)
        gsInfo << "Least squares finished" << std::endl;
    // --- Smoothing matrix computation
    assembleContinuous(entries_mat, m_B);

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    if(m_print_messages)
        gsInfo << "Assembling finished" << std::endl;

    m_A_mat.resize(num_basis, num_basis);
    m_A_mat.setZero();
/*
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_basis->dim(); ++i) // to do: improve
        //nonZerosPerCol *= 2 * m_basis->degree(i) + 1;
        nonZerosPerCol *= ( 2 * m_basis->degree(i) + 1 ) * 4;//AM: Related to bandwidth of basis
    m_A_mat.reserve( gsVector<index_t>::Constant(num_basis, nonZerosPerCol ) );
*/
    m_A_mat.setFrom(entries_mat);
    m_A_mat.makeCompressed();

    typename gsSparseSolver<T>::BiCGSTABILUT solver( m_A_mat );

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
        m_result_linear = NULL;
        return;
    }

    // Solves for many right hand side  columns
    gsMatrix<T> x, x2;
    gsMatrix<T>* x_ref = &x;
    x = solver.solve(m_B);
    if(m_print_messages)
        gsInfo << "Inverse computed" << std::endl;

    /// In case we could not split the dimensions for the resolution,
    /// we now split the dimensions in the vector x
    if(! can_split)
    {
        num_basis = sizeDof();
        x2 = gsMatrix<T>(num_basis, m_dim);
        for(int i = 0;i < num_basis;i++)
        {
            for(unsigned d1 = 0;d1 < m_dim;d1++)
            {
                x2(i, d1) = x(i*m_dim + d1);
            }
        }
        x_ref = &x2;
    }
    GISMO_ASSERT(m_result_linear == NULL, "The result should be NULL before the computation");
    _this->makeGeometry(*x_ref);
}


template<short_t d, class GenBasis, class T>
bool _gsFittingLinSyst<d, GenBasis, T>::has_smoothing()
{
    if(m_quadr_1d != NULL && m_quadr_1d->has_smoothing())
        return true;
    if(m_quadr_2d != NULL && m_quadr_2d->has_smoothing())
        return true;
    if(m_quadr_3d != NULL && m_quadr_3d->has_smoothing())
        return true;
    return false;
}


template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
assembleContinuous(gsSparseEntries<T> & entries_mat,
                   gsMatrix<T>& m_B)
{
    if(m_quadr_1d != NULL)
        m_quadr_1d->assembleSystem(entries_mat, m_B,
                                   m_split_dim, m_isLinear);
    if(m_quadr_2d != NULL)
        m_quadr_2d->assembleSystem(entries_mat, m_B,
                                   m_split_dim, m_isLinear);
    if(m_quadr_3d != NULL)
        m_quadr_3d->assembleSystem(entries_mat, m_B,
                                   m_split_dim, m_isLinear);
}

template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
set_current_map(gsFunctionSet<T>* current)
{
    if(m_current != NULL && m_current != m_template)
        delete m_current;
    m_current = current;
    if(m_least_squares != NULL)
        m_least_squares->set_current_map(current);
    if(m_quadr_1d != NULL)
        m_quadr_1d->set_current_map(current);
    if(m_quadr_2d != NULL)
        m_quadr_2d->set_current_map(current);
    if(m_quadr_3d != NULL)
        m_quadr_3d->set_current_map(current);
}

template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
decrease_coeff_smoothing(T ratio)
{
    if(m_quadr_1d != NULL)
        m_quadr_1d->decrease_coeff_smoothing(ratio);
    if(m_quadr_2d != NULL)
        m_quadr_2d->decrease_coeff_smoothing(ratio);
    if(m_quadr_3d != NULL)
        m_quadr_3d->decrease_coeff_smoothing(ratio);
}

template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
decrease_linear_coeff_smoothing(T ratio)
{
    if(m_quadr_1d != NULL)
        m_quadr_1d->decrease_linear_coeff_smoothing(ratio);
    if(m_quadr_2d != NULL)
        m_quadr_2d->decrease_linear_coeff_smoothing(ratio);
    if(m_quadr_3d != NULL)
        m_quadr_3d->decrease_linear_coeff_smoothing(ratio);
}


template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::makeGeometry(gsMatrix<T>& coefs)
{
    m_result_linear = makeGeometryGen<T>
        (*m_basis, coefs, m_local_global);
}

/*
template<short_t d, class T>
void gsFittingLinSyst<d, gsMultiBasis<T>, T>::
makeGeometry(gsMatrix<T>& coefs)
{
    m_result_linear = makeGeometryGen<T>
        (m_basis, coefs, m_local_global);
}
*/


template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
addGlobalIntegrand(gsFittingIntegrand<d,T>* integr)
{
    gsFittingLinSyst<d, GenBasis, T>* _this = static_cast
        <gsFittingLinSyst<d, GenBasis, T>*>(this);
    addGlobalIntegrandGen<GenBasis, T>(_this, integr);
}

template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
add_regularization_term(T coeff)
{
    this->addGlobalIntegrand
        (new gsFittingIntegrandL2Dist<d, gsGeometry<T>, T>
         (m_dim, NULL, false, coeff) );
    this->addGlobalIntegrand
        (new gsFittingIntegrandLin<d, 2, T>
         (coeff, 1., 0., false) );
}


template<short_t d, class GenBasis, class T>
 int _gsFittingLinSyst<d, GenBasis, T>::
add_dimension_to_index(int global_index)
{
    if(m_split_dim)
        return global_index;
    else
        return d * global_index;
}


template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
actualize_basis()
{

    gsFittingLinSyst<d, GenBasis, T>* _this = static_cast
        <gsFittingLinSyst<d, GenBasis, T>*>(this);

    _this->actualize_multipatch(true);
    if(m_least_squares != NULL)
        m_least_squares->actualize_degrees_freedomLS();
    if(m_quadr_1d != NULL)
        m_quadr_1d->actualize_degrees_freedom();
    if(m_quadr_2d != NULL)
        m_quadr_2d->actualize_degrees_freedom();
    if(m_quadr_3d != NULL)
        m_quadr_3d->actualize_degrees_freedom();
}

template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
actualize_identity_mapping()
{
    if(m_quadr_1d != NULL)
        m_quadr_1d->set_identity_mapping(m_template);
    if(m_quadr_2d != NULL)
        m_quadr_2d->set_identity_mapping(m_template);
    if(m_quadr_3d != NULL)
        m_quadr_3d->set_identity_mapping(m_template);
}

/// <------------- Energy computation --------------->


template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
print_error(bool compute)
{
    if(compute)
        computeErrors();
    m_energy_max = maxEnergySmoothing();
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
            gsInfo << "ENERGY : max : " << maxEnergySmoothing()
                   << "   min : " << minEnergySmoothing()
                   << "    total : " << totalEnergySmoothing()
                   << std::endl;
        }
    }
}

template<short_t d, class GenBasis, class T>
T _gsFittingLinSyst<d, GenBasis, T>::
maxError()
{
    T res = 0;
    if(m_least_squares != NULL)
        res = m_least_squares->maxError();
    if(m_quadr_1d != NULL)
        res = std::max(res, m_quadr_1d->maxError());
    if(m_quadr_2d != NULL)
        res = std::max(res, m_quadr_2d->maxError());
    if(m_quadr_3d != NULL)
        res = std::max(res, m_quadr_3d->maxError());
    return res;
}

template<short_t d, class GenBasis, class T>
T _gsFittingLinSyst<d, GenBasis, T>::
minError()
{
    T res = std::numeric_limits<T>::max();
    if(m_least_squares != NULL)
        res = m_least_squares->minError();
    if(m_quadr_1d != NULL)
        res = std::min(res, m_quadr_1d->minError());
    if(m_quadr_2d != NULL)
        res = std::min(res, m_quadr_2d->minError());
    if(m_quadr_3d != NULL)
        res = std::min(res, m_quadr_3d->minError());
    return res;
}

template<short_t d, class GenBasis, class T>
T _gsFittingLinSyst<d, GenBasis, T>::
totalError()
{
    T res = 0;
    if(m_least_squares != NULL)
        res = m_least_squares->totalError();
    if(m_quadr_1d != NULL)
        res += m_quadr_1d->totalError();
    if(m_quadr_2d != NULL)
        res += m_quadr_2d->totalError();
    if(m_quadr_3d != NULL)
        res += m_quadr_3d->totalError();
    return res;
}


template<short_t d, class GenBasis, class T>
T _gsFittingLinSyst<d, GenBasis, T>::
maxEnergySmoothing()
{
    T res = 0;
    T val;
    if(m_quadr_1d != NULL)
    {
        val = m_quadr_1d->maxEnergySmoothing();
        if(val == -1.)
            return val;
        res += val;
    }
    if(m_quadr_2d != NULL)
    {
        val = m_quadr_2d->maxEnergySmoothing();
        if(val == -1.)
            return val;
        res += val;
    }
    if(m_quadr_3d != NULL)
    {
        val = m_quadr_3d->maxEnergySmoothing();
        if(val == -1.)
            return val;
        res += val;
    }
    return res;
}

template<short_t d, class GenBasis, class T>
T _gsFittingLinSyst<d, GenBasis, T>::
minEnergySmoothing()
{
    T res = 0.;
    T val;
    if(m_quadr_1d != NULL)
    {
        val = m_quadr_1d->minEnergySmoothing();
        if(val == -1.)
            return val;
        res += val;
    }
    if(m_quadr_2d != NULL)
    {
        val = m_quadr_2d->minEnergySmoothing();
        if(val == -1.)
            return val;
        res += val;
    }
    if(m_quadr_3d != NULL)
    {
        val = m_quadr_3d->minEnergySmoothing();
        if(val == -1.)
            return val;
        res += val;
    }
    return res;
}

template<short_t d, class GenBasis, class T>
T _gsFittingLinSyst<d, GenBasis, T>::
totalEnergySmoothing()
{
    T res = 0;
    T val;
    if(m_quadr_1d != NULL)
    {
        val = m_quadr_1d->totalEnergySmoothing();
        if(val == -1.)
            return val;
        res += val;
    }
    if(m_quadr_2d != NULL)
    {
        val = m_quadr_2d->totalEnergySmoothing();
        if(val == -1.)
            return val;
        res += val;
    }
    if(m_quadr_3d != NULL)
    {
        val = m_quadr_3d->totalEnergySmoothing();
        if(val == -1.)
            return val;
        res += val;
    }
    return res;
}

template<short_t d, class GenBasis, class T>
void _gsFittingLinSyst<d, GenBasis, T>::
computeErrors(int type)
{
    GISMO_ASSERT(m_current != NULL,
                 "Error: there is no solution");
    if(m_least_squares != NULL)
        m_least_squares->computeErrors(type);
    if(m_quadr_1d != NULL)
        m_quadr_1d->computeEnergies(type);
    if(m_quadr_2d != NULL)
        m_quadr_2d->computeEnergies(type);
    if(m_quadr_3d != NULL)
        m_quadr_3d->computeEnergies(type);
}

template<short_t d, class GenBasis, class T>
gsFittingLinSyst<d, GenBasis, T>::
gsFittingLinSyst(unsigned dim_im, GenBasis & basis,
                 bool remove_result,
                 gsFittingParam<T>& param,
                 gsLeastSquares<GenBasis, T>*
                 least_squares) :
Base::_gsFittingLinSyst(dim_im, basis, remove_result,
                        param, least_squares){ }


template<short_t d, class T>
gsFittingLinSyst<d, gsBasis<T>, T>::
gsFittingLinSyst(unsigned dim_im, gsBasis<T> & basis,
                 bool remove_result,
                 gsFittingParam<T>& param,
                 gsLeastSquares<gsBasis<T>, T>*
                 least_squares) :
Base::_gsFittingLinSyst(dim_im, basis, remove_result,
                        param, least_squares)  { }


template<short_t d, class T>
gsFittingLinSyst<d, gsMultiBasis<T>, T>::
gsFittingLinSyst(unsigned dim_im,
                 gsMultiBasis<T> & basis,
                 bool remove_result,
                 gsFittingParam<T>& param,
                 gsLeastSquares<gsMultiBasis<T>, T>*
                 least_squares) :
Base::_gsFittingLinSyst(dim_im, basis, remove_result,
                        param, least_squares){  }

} // namespace gismo

/** @file gsFittingIterative.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/

#pragma once


#include <gsFitting/gsFittingQuadrature.hpp>
#include <gsFitting/gsFittingIdConstr.hpp>
#include <gsFitting/gsLeastSquares.hpp>
#include <gsFitting/gsFittingIterative.h>
#include <gsFitting/gsFittingUtilsProjection.h>
#include <gsFitting/gsFittingUtilsIO.hpp>
#include <gsFitting/gsPointsSampling.hpp>

#define ERROR_DETER_NEG gsWarn << "The determinant is somewhere negative!! It seems like the basis should be refined or that a higher error should be expected. The constraints imposed seem too strong." << std::endl

#define ERROR_DETER_NEG_PROJ gsWarn << "Problem during the projection! The determinant is somewhere negative after projection!! It seems like the basis should be refined or that a higher error should be expected. The constraints imposed seem too strong." << std::endl

namespace gismo
{

template<short_t d, class GenBasis, class T>
void gsFittingIter<d, GenBasis, T>::compute()
{
    /***** Before the nonlinear loop ****/
    initComputation();

    if(m_print_messages)
        gsInfo << "----------COMPUTE RESULT----------" << std::endl;

    bool finished = false;
    bool reinitialized, cancelled;
    T last_ener_max, energyMax;

    int frequ_init = m_frequ_proj / 2;

    m_ind_iter = 0;
    m_global_ind = 0;
    while(! finished)
    {
        reinitialized = false;
        if(m_print_messages)
            gsInfo << "-------------- ITERATION "
                   << util::to_string(m_ind_iter) << " ---------------"
                   << std::endl;
        this->computeLinear();
        GISMO_ASSERT(m_result_linear != NULL,
                     "Error during the construction of the deformation");

        if(m_is_displacement)
        {
            computeNormDisplacement();
            if(m_print_messages)
                gsInfo << "norm displacement : L2: "
                       << m_displacement << "   L^inf: "
                       << m_displacement_Linf << std::endl;
        }
        /*
          if(!m_isLinear && m_afterLinearSolving != NULL)
          reinitialized = m_afterLinearSolving
          ->iterAdaptFitting(*this, m_ind_iter); */

        /// reinitialized is always true
        if(! reinitialized)
        {
            cancelled = false;
            if(m_is_displacement)
            {
                localResulToGlobalResult();
                last_ener_max = this->maxEnergySmoothing();
                this->computeErrors();

                energyMax = this->maxEnergySmoothing();
                if(energyMax < 0)
                {
                    gsWarn << "Negative energy" << std::endl;
                    cancelled = true;
                }
                else
                    cancelled = ! validDisplacement
                        (energyMax, last_ener_max);
                if(cancelled)
                    finished = cancel_displacement();
                if(energyMax < 0 && finished)
                    ERROR_DETER_NEG;
            }
            else
            {
                GISMO_ASSERT(this->isLinear(),
                             "The system is nonlinear. It should be linearized");
                this->set_current_map(m_result_linear);
            }
            if(! cancelled)
            {
                if(m_use_lagrangeM)
                {
                    if(m_ind_iter % m_frequ_proj == frequ_init)
                    {
                        this->print_error(! m_is_displacement);
                        exportIteration("_beforeProj");
                        cancelled = ! applyLagrangeMMethod();
                        if(cancelled)
                        {
                            if(m_iter != NULL)
                            {
                                finished = m_iter->
                                    after_cancelled_iter
                                    (*this, m_ind_iter);
                                m_ind_iter++;
                            }
                            else
                                finished = true;
                            if(finished)
                                ERROR_DETER_NEG_PROJ;
                        }
                    }
                }
                if(! cancelled)
                {
                    this->print_error();
                    if(m_print_messages)
                        gsInfo << std::endl;

                    exportIteration();

                    finished = isFinished(m_ind_iter);

                    if(finished)
                    {
                        printEndIterations();
                        if(m_iter != NULL)
                            finished = m_iter->
                                iterAdaptFitting(*this, m_ind_iter);

                        if(! finished)
                        {
                            m_ind_iter = 0;
                            m_global_ind++;
                        }
                    }
                    else
                        m_ind_iter++;

                    if(m_print_messages)
                        gsInfo << std::endl << std::endl << std::endl;
                }
            }
        }
        if(! finished || m_is_displacement)
        {
            if(! m_is_displacement)
                this->set_current_map(NULL);
            else
                delete m_result_linear;
            m_result_linear = NULL;
        }
    }
    if(m_print_messages)
        gsInfo << "------------ END PROCESS --------------"
               << std::endl;

}


template<short_t d, class GenBasis, class T>
void gsFittingIter<d, GenBasis, T>::
computeNormDisplacement()
{
    gsFittingQuadrature<d, GenBasis, 2, T>
        L2Ener(m_basis, m_dim);
    L2Ener.addIntegrand(new gsFittingIntegrandL2Norm<d, T>() );
    L2Ener.set_current_map(m_result_linear);

    L2Ener.computeEnergies();
    m_displacement = L2Ener.totalError();
    m_displacement_Linf = L2Ener.maxError();
}



template<short_t d, class GenBasis, class T>
bool gsFittingIter<d, GenBasis, T>::isFinished(int iter)
{
    if(iter >= m_Nmax_iter - 1)
        return true;
    /// If the iterative method does not compute a
    /// displacement at each step.
    /// The stopping criterion must be implemented. Otherwise we suppose
    /// that the system is linear and we only do 1 iteration.
    if(! m_is_displacement)
        return true;
    return m_displacement < m_threshold_displacement;
}


template<short_t d, class GenBasis, class T>
void gsFittingIter<d, GenBasis, T>::
localResulToGlobalResult()
{
    GISMO_ASSERT(m_current != NULL && m_result_linear != NULL,
                 "Error: either the current result is NULL or the local linear solution is NULL");
    sumGenericGeometries(m_basis, *m_current,
                         *m_result_linear);
}



template<short_t d, class GenBasis, class T>
gsFunctionSet<T>* gsFittingIter<d, GenBasis, T>::
refineIdentityMapping(bool sampling)
{
    if(sampling)
        return basisProjectionSampling<d, GenBasis, T>
            (m_basis, m_dim, *m_identity);
    else
        return basisProjectionContinuous<d, GenBasis, T>
            (m_basis, m_dim, *m_identity);
}


template<short_t d, class GenBasis, class T>
void gsFittingIter<d, GenBasis, T>::
initComputation()
{
    GISMO_ASSERT(m_current == NULL,
                 "Result already computed");

    this->set_properties_integrands();

    if(!this->isLinear() || m_use_lagrangeM)
        m_is_displacement = true;
    else if(! m_use_lagrangeM)
        m_is_displacement = false;

    /// Construct the "identity mapping"
    bool use_template = ( m_is_displacement && m_initial_map == NULL)
        || m_deform_minim;
    if(use_template && m_identity == NULL)
    {
        this->constructId();
        GISMO_ASSERT(m_identity != NULL, "Error during the construction of the template");
        if(m_deform_minim)
            this->actualize_identity_mapping();
    }

    if(m_initial_map == NULL && m_is_displacement)
    {
        m_initial_map = refineIdentityMapping();
        GISMO_ASSERT(m_initial_map != NULL, "Error during the initialization of the algorithm");
    }
    if(m_is_displacement)
        this->set_current_map(m_initial_map);

    if(m_export_initialization && m_is_displacement)
    {
        std::string parameter_lines_filename
            = m_output + "_template";
        writeParameterLines(getGeometry(m_basis, *m_initial_map),
                            parameter_lines_filename,
                            50, m_num_lines);
    }

}

template<short_t d, class GenBasis, class T>
void gsFittingIter<d, GenBasis, T>::exportNodes()
{
    if(m_export_nodes)
    {
        unsigned s_patch = nBasesGen(m_basis);

        for(unsigned patch = 0;patch < s_patch;patch++)
        {
            std::string parameter_lines_filename
                = m_output + "_globalInd_"
                + util::to_string(m_global_ind)
                + "_patch" + util::to_string(patch);
            writeKnots(getNthGeometry(m_basis, *m_current,
                                      patch),
                       parameter_lines_filename);
        }
    }
}

template<short_t d, class GenBasis, class T>
bool gsFittingIter<d, GenBasis, T>::
applyLagrangeMMethod()
{
    if(m_use_lagrangeM)
    {
        gsFunctionSet<T>* func_copy = m_current->clone().release();
        gsFunctionSet<T>* proj = projectSolution();
        this->computeErrors();
        if(this->maxEnergySmoothing() < 0)
        {
            this->set_current_map(func_copy);
            return false;
        }
        GISMO_ASSERT(proj == m_current,
                     "The two mappings should be the same");
        m_energy.adaptMultipliers(func_copy, proj/*, m_A_mat*/);
        return true;
    }
    return true;
}

template<short_t d, class GenBasis, class T>
gsFunctionSet<T>* gsFittingIter<d, GenBasis, T>::
projectSolution()
{
    gsFittingParam<T> param;
    bool remove_result = false;
    T reg = 50. * math::pow(m_displacement, 2);
    gsInfo << "Tikhonov coeff: " << reg << std::endl;

    /// We set the parameters of the computation
    param.coeff_NL_global = 0.;
    param.use_lagrangeM = false;
    param.export_iterations = false;
    param.max_num_iter_NL = 2;
    param.p_identity = m_identity;
    param.p_initialization = m_current;
    param.export_points = false;
    param.output = m_output + "_iter_" +
        util::to_string(m_ind_iter) + "_proj";
    param.threshold_NL = 0.00000001;

    param.coeff_linear_global = reg;
    param.coeff_linear_gradient = 1.;
    param.coeff_linear_hessian = 0.;
    param.deform_min = false;
    param.coeff_NL_global = 0.;

    gsFittingIterNL <d, GenBasis, T>
        fitting(m_dim, m_basis, remove_result, param);
    fitting.energy().resetPointsLS
        (this->leastSquares());

    fitting.add_regularization_term(reg);
    fitting.compute();

    return fitting.result();
}

template<short_t d, class GenBasis, class T>
void gsFittingIter<d, GenBasis, T>::
exportIteration(const std::string& add_name)
{
    if(m_export_iterations)
    {
        int nPoints = 50;
        GISMO_ASSERT(m_current != NULL,
                     "Error current solution NULL");
        std::string filename
            = m_output + add_name + "_globalInd_"
            + util::to_string(m_global_ind)
            + "_iter_" + util::to_string(m_ind_iter);
        std::string parameter_lines_filename
            = filename + "_full";
        writeParameterLines(getGeometry(m_basis, *m_current),
                            parameter_lines_filename,
                            nPoints, m_num_lines);

        if(d == 3 && m_nSlices > 0)
        {
            int dir_slices = 2;
            exportSlices(getGeometry(m_basis, *m_current),
                         filename, dir_slices, m_nSlices, nPoints);
        }
    }
}

template<short_t d, class GenBasis, class T>
void gsFittingIter<d, GenBasis, T>::printEndIterations()
{
    if(m_print_messages)
    {
        gsInfo << std::endl << std::endl <<
            "-----------------------------------------------------"
               << std::endl;
        gsInfo
            << "-------- END OF NONLINEAR ITERATIONS NUM "
            << m_global_ind << "---------" << std::endl;
        gsInfo <<
            "-----------------------------------------------------"
               << std::endl << std::endl << std::endl;
    }
}


template<short_t d, class GenBasis, class T>
void gsFittingIter<d, GenBasis, T>::
reinitialize(gsFunctionSet<T>* new_init)
{
    m_initial_map = new_init;
    if(m_current != m_identity && m_current != NULL)
        delete m_current;
    this->set_current_map(m_initial_map);
}

template<short_t d, class GenBasis, class T>
bool gsFittingIter<d, GenBasis, T>::
cancel_displacement()
{
    if(m_print_messages)
        gsInfo << "We cancel the last displacement"
               << std::endl;
    sumGenericGeometries(m_basis, *m_current,
                         *m_result_linear, -1.);
    if(m_ind_iter != 0)
    {
        if(m_iter != NULL)
            return m_iter->after_cancelled_iter(*this, m_ind_iter);
        /// The displacement is too important
        /*   gsWarn << "The displacement is too important. This case should be considered (adding Tikhonov regularization?). Most probably the algorithm will not converge."
             << std::endl;*/
    }
    else
    {
        /// We just modified the optimization problem.
        /// This modification must be cancelled
        GISMO_ASSERT(m_iter != NULL, "Should not be here. How the problem can be modified if m_iter is NULL?");
        return m_iter->cancel(*this, m_ind_iter);
    }
    return true;
}

template<short_t d, class GenBasis, class T>
bool gsFittingIter<d, GenBasis, T>::
validDisplacement(T newEner, T lastEner)
{
    return true;
    if(m_ind_iter == 0 && m_global_ind == 0)
        return true;
    if(newEner > 1.33 * lastEner)
    {
        if(m_print_messages)
            gsInfo << "last energy: " << lastEner
                   << "   new energy: " << newEner << std::endl;
        return false;
    }
    else
        return true;
}

} // namespace gismo

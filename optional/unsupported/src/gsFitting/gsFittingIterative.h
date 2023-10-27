/** @file gsFittingIter.h

    @brief Contains the class gsFittingIter providing iterative
    methods for splitting. Used for hierarchical fitting, fitting with
    nonlinear smoothing energies,...

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <gsFitting/gsParameterLines.h>
#include <gsFitting/gsFittingUtilsIO.h>

#include <gsMatrix/gsMatrix.h>
#include <gsFitting/gsPointsSampling.h>

#include <gsFitting/gsFittingIdConstr.h>
#include <gsFitting/gsFittingQuadrature.h>
#include <gsFitting/gsLeastSquares.h>

namespace gismo
{

/**
   @brief Interface used for iterative fitting. Contains some function
   that can be called after solving the linear system at each iteration.
   Can be used for example for
   hierarchical refinement or geometry fitting with adaptative parametrization
*/
template <short_t d, class GenBasis, class T>
class gsFittingImproveAfterCV;


/**
   @brief
   Class for performing an iterative algorithm for
   fitting some data

   GenBasis the type of basis used (only gsBasis<T> and gsMultiBasis<T> are used)
**/
template<short_t d, class GenBasis, class T=real_t>
class gsFittingIter : public gsFittingBase<d, GenBasis, T>
{
    public:
    typedef gsFittingBase<d, GenBasis, T> Base;

    /// Constructors
    explicit gsFittingIter(bool is_displacement, int Nmax_iter,
                           T threshold_displacement,
                           unsigned dim_im, GenBasis& basis,
                           bool remove_result,
                           gsFittingParam<T>& param,
                           gsFittingImproveAfterCV<d, GenBasis, T>*
                           iter = NULL,
                           gsFittingImproveAfterCV<d, GenBasis, T>*
                           afterLinearSolving = NULL) :
    Base::gsFittingBase(dim_im, basis, remove_result, param)
    {
        init_iter(is_displacement, Nmax_iter,
                  threshold_displacement, param, iter,
                  afterLinearSolving);
    }

    explicit gsFittingIter(bool is_displacement, int Nmax_iter,
                           T threshold_displacement,
                           unsigned dim_im, GenBasis& basis,
                           bool remove_result,
                           gsFittingParam<T>& param,
                           gsBoundaryConditions<T>& bc,
                           gsLeastSquares<GenBasis, T>*
                           least_squares = NULL,
                           gsFittingImproveAfterCV<d, GenBasis, T>*
                           iter = NULL,
                           gsFittingImproveAfterCV<d, GenBasis, T>*
                           afterLinearSolving = NULL) :
    Base::gsFittingBase(dim_im, basis, remove_result,
                           param, bc, least_squares)
    {
        init_iter(is_displacement, Nmax_iter,
                  threshold_displacement, param, iter,
                  afterLinearSolving);
    }

    void init_iter(bool is_displacement, int Nmax_iter,
                   T threshold_displacement,
                   gsFittingParam<T>& param,
                   gsFittingImproveAfterCV<d, GenBasis, T>*
                   iter = NULL,
                   gsFittingImproveAfterCV<d, GenBasis, T>*
                   afterLinearSolving = NULL)
    {
        m_is_displacement = is_displacement;
        m_initial_map = param.p_initialization;
        m_threshold_displacement = threshold_displacement;
        m_Nmax_iter = Nmax_iter;
        m_iter = iter;
        m_displacement = 0.;
        m_displacement_Linf = 0.;
        m_displacement_max = 0.;
        m_afterLinearSolving = afterLinearSolving;
        m_use_lagrangeM = param.use_lagrangeM;

        m_export_nodes = param.export_nodes;
        m_export_initialization = param.export_initialization;

        m_num_lines = param.num_lines;

        m_nSlices = param.nSlices;

        m_ind_iter = 0;
        m_global_ind = 0;

        m_export_iterations = param.export_iterations;

        m_frequ_proj = param.frequ_proj;
    }
    /// Adds the solution of the linearized system to
    /// the current solution
    void localResulToGlobalResult();

    /// Returns true if the iterative algorithm is finished.
    /// Here, the algorithm is finished if the maximal number
    /// of iterations have been reached or if the norm of the
    /// displacement is below some threshold
    virtual bool isFinished(int iter);

    /// Performs the iterative fitting
    /// The result is set in m_current
    void compute();

    /// Computes the norm of the current displacement.
    /// the result is set in m_displacement
    void computeNormDisplacement();

    /// Returns the computed fitting
    gsFunctionSet<T>* result() const
    { return m_current; }

    /// Returns true if the algorithm consists in displacing
    /// the current solution until some criterion is reached.
    bool is_displacement()
    {  return m_is_displacement;  }

    /// Returns the norm of the displacement
    T normDisplacement()
    {  return m_displacement;  }

    /// Returns the geometry of the identity mapping defined on
    /// the basis m_basis
    gsFunctionSet<T>* refineIdentityMapping(bool sampling = false);

    /// Returns the geometry used for the initialization
    /// of the algorithm
    gsFunctionSet<T>* initialPos()
    {   return m_initial_map;   }

    /// Prints some information. Called in the end of each iteration
    void printEndIterations();

    /// Export the nodes of the basis and its image
    /// by the current solution
    void exportNodes();

    /// Sets the initial map to new_init.
    void reinitialize(gsFunctionSet<T>* new_init);

    /// Exports the result. Called in the end of each iteration
    void exportIteration(const std::string& add_name = "");

    /// Initialize some data (computes the identity mapping if
    /// needed,...) before starting the iterations.
    void initComputation();

    /// Project the solution on the constraint (the data to be fitted).
    /// This projection is performed by removing the smoothing energy
    /// and adding a Tikhonov regularization.
    /// Only implemented for the moment for Least Squares method.
    gsFunctionSet<T>* projectSolution();

    /// Called after having computed the solution linearized energy.
    /// Project the new solution and adapt the Lagrange multipliers
    /// Return false if this step generated a mapping that is not inversible
    bool applyLagrangeMMethod();

    /// In this case the displacement is not valid, we cancel it and
    /// modify the energy so that the next displacement is smaller
    bool cancel_displacement();

    /// Can the displacement be considered as valid (not to big)
    bool validDisplacement(T newEner, T lastEner);

public:
    /// Pointer to the initial mapping (geometry)
    gsFunctionSet<T>* m_initial_map;

    /// We stop when the norm of the displacement is below this threshold
    /// Used only if we compute a displacement at each step
    T m_threshold_displacement;

    /// The maximal number of iterations
    int m_Nmax_iter;

    /// Performs the different computations needed in the end of each iteration.
    /// Interface that can be used for example for hierarchical refinement (done) or
    /// modification of the parametrization of the target geometry
    gsFittingImproveAfterCV<d, GenBasis, T>* m_iter;

    /// Called after the computation of the solution of the linear system.
    /// In the current implementation, it is only used to
    /// recompute the initialization in the case where
    /// the displacement is too important in the first step
    /// (NOT USED NOW)
    gsFittingImproveAfterCV<d, GenBasis, T>* m_afterLinearSolving;

    /// The norm (L2) of the displacement made in the last step
    T m_displacement;
    /// The norm (L^inf) of the displacement made in the last step
    T m_displacement_Linf;

    /// The norm of the displacement made in the last step
    T m_displacement_max;

    /// Do we export the nodes after each reffinement?
    bool m_export_nodes;

    /// Do we export the geometry used for the initialization of the algorithm?
    bool m_export_initialization;

    /// Apply the algorithm that consists in computing Lagrange
    /// multipliers if true
    bool m_use_lagrangeM;

    /// The index of the current iteration
    int m_ind_iter;

    /// The global index of the current iteration
    int m_global_ind;

    /// The number of slices used for the exportation in case d = 3
    int m_nSlices;

    /// The number of coordinate curves used for the exportation
    int m_num_lines;

    /// True if each iteration must be exported
    bool m_export_iterations;

    /// If the projection algorithm is applied,
    /// frequency of its application. If one, it is applied
    /// every single iteration. If two, it is applied one over
    /// two iterations...
    int m_frequ_proj;

    using Base::m_current;
    using Base::m_result_linear;
    using Base::m_identity;
    using Base::m_print_messages;
    using Base::m_output;
    using Base::m_basis;
    using Base::m_deform_minim;
    using Base::m_is_displacement;
    using Base::m_energy;

    using Base::m_dim;

};  /// class gsFittingIter



/**
   Class used for nonlinear fitting
**/
template<short_t d, class GenBasis, class T=real_t>
class gsFittingIterNL : public gsFittingIter<d, GenBasis, T>
{
    public:
    typedef gsFittingBase<d, GenBasis, T> Base;

    /// Constructors
    explicit gsFittingIterNL(unsigned dim_im, GenBasis& basis,
                             bool remove_result,
                             gsFittingParam<T>& param,
                             gsFittingImproveAfterCV<d, GenBasis, T>*
                             iter = NULL,
                             gsFittingImproveAfterCV<d, GenBasis, T>*
                             afterLinearSolving = NULL) :
    gsFittingIter<d, GenBasis, T>::
    gsFittingIter(true, param.max_num_iter_NL, param.threshold_NL,
                  dim_im, basis, remove_result, param,
                  iter, afterLinearSolving)
    { }


public:
};  /// class gsFittingIterNL

/**
   @brief Interface used for iterative fitting. Contains some function
   that can be called udring the iterative process. Only called for the moment
   once the process have converged for hierarchical refinement and to decrease
   the smoothing coefficients when necessary.
   Could also be used for adapting the parametrization of the boundary.
*/
template <short_t d, class GenBasis, class T>
class gsFittingImproveAfterCV
{
public:
    gsFittingImproveAfterCV(){ }

    /// Performs operations at the end of the nonlinear loop
    /// Returns false if we should restart the nonlinear process
    virtual bool iterAdaptFitting(gsFittingIter<d, GenBasis, T>&
                                  fitting, int iter) = 0;

    /// In case the displacement is too important after calling
    /// iterAdaptFitting. This function is called.
    virtual bool cancel(gsFittingIter<d, GenBasis, T>&
                        fitting, int ind_iter)
    { GISMO_NO_IMPLEMENTATION; }

    /// In case one iteration is too important,
    /// this function is called.
    virtual bool after_cancelled_iter(gsFittingIter<d, GenBasis, T>&
                                      fitting, int ind_iter)
    { GISMO_NO_IMPLEMENTATION; }
};

/*
 Class used for decreasing the coefficients of the smoothing energy once
 the system has converged.
 */
template <short_t d, class GenBasis, class T>
class gsFittingAdaptSmoo :
        public gsFittingImproveAfterCV<d, GenBasis, T>
{
public:
    gsFittingAdaptSmoo(T coeff, T coeff_min, int ind_max,
                       T coeff_init_lin, T coeff_init_NL,
                       bool print_messages,
                       gsFittingImproveAfterCV<d, GenBasis, T>*
                       other_call = NULL)
    {
        m_ind = 0;
        m_ind_max = ind_max;
        m_coeff_reduc = coeff;
        m_coeff_reduc_init = coeff;
        m_coeff_min = coeff_min;
        m_other_call = other_call;
        m_print_messages = print_messages;
        m_coeff_tot = 1.;
        m_coeff_tot_lin = 1.;
        m_coeff_init_lin = coeff_init_lin;
        m_coeff_init_NL = coeff_init_NL;

        m_ratio = 1.;
    }

    /// Performs operations at the end of the nonlinear loop
    /// Returns false if we should restart the nonlinear process
    bool iterAdaptFitting(gsFittingIter<d, GenBasis, T>&
                          fitting, int iter)
    {
        /// In case there have only been one iteration
        /// and all the conditions are satisfied,
        /// we decrease the smoothing coefficient
        gsInfo << "new coeff lin: " << m_coeff_tot_lin * m_coeff_init_lin * m_coeff_reduc;
        gsInfo << "   new coeff NL: " << m_coeff_tot * m_coeff_init_NL * m_coeff_reduc
               << "    coeff min: " << m_coeff_min << std::endl;
        gsInfo << "coeff: " << m_coeff_reduc << "   ind: " << m_ind
               << "  max: " << m_ind_max << std::endl;

        if(iter > 0 && m_ind < m_ind_max
           && ( m_coeff_tot * m_coeff_init_NL * m_coeff_reduc > m_coeff_min
                || m_coeff_init_lin * m_coeff_tot_lin * m_coeff_reduc
                > m_coeff_min) )
        {
            if(m_coeff_init_lin * m_coeff_tot_lin
               > m_coeff_init_NL * m_ratio){
                fitting.decrease_linear_coeff_smoothing(m_coeff_reduc);
                m_coeff_tot_lin *= m_coeff_reduc;
            }
            else{
                fitting.decrease_coeff_smoothing(m_coeff_reduc);
                m_coeff_tot *= m_coeff_reduc;
                m_coeff_tot_lin *= m_coeff_reduc;
            }
            fitting.decreaseLMult(m_coeff_reduc);
            m_ind++;
            //         m_coeff_reduc /= 1.5;
            if(m_print_messages)
                gsInfo << "The coefficient has been multiplied by: "
                       << m_coeff_tot << std::endl;
            return false;
        }
        if(m_other_call != NULL)
        {
            m_ind = 0.;
            return m_other_call->
                iterAdaptFitting(fitting, iter);
        }
        else return true;
    }

    bool cancel(gsFittingIter<d, GenBasis, T>& fitting,
                int ind_iter)
    {
        if(m_coeff_init_lin * m_coeff_tot_lin
           > m_coeff_init_NL * m_ratio){
            fitting.decrease_linear_coeff_smoothing(1./m_coeff_reduc);
            m_coeff_tot_lin /= m_coeff_reduc;
        }
        else{
            fitting.decrease_coeff_smoothing(1./m_coeff_reduc);
            m_coeff_tot /= m_coeff_reduc;
        }
        m_coeff_reduc = (1. + m_coeff_reduc)/2.;
        if(m_coeff_reduc > 0.9)
        {
            m_coeff_reduc = m_coeff_reduc_init;
            return after_cancelled_iter(fitting, ind_iter);
        }
        iterAdaptFitting(fitting, ind_iter+1);
        return false;
    }


    bool after_cancelled_iter(gsFittingIter<d, GenBasis, T>&
                              fitting, int ind_iter)
    {
        if(m_other_call != NULL)
            return m_other_call->iterAdaptFitting(fitting, ind_iter);
        return true;
    }

protected:
    /// The number of time iterAdaptFitting have already been called.
    int m_ind;

    /// The maximal number of time the coefficients can be decreased
    int m_ind_max;

    /// The minimal smoothing coefficient that can be set.
    T m_coeff_min;

    /// We multiply the smoothing coefficients by this coefficient
    /// once the process have converged.
    T m_coeff_reduc;

    /// The coefficient initially set.
    T m_coeff_reduc_init;

    /// The total reduction operated
    T m_coeff_tot;

    /// The total reduction operated on linear coefficients
    T m_coeff_tot_lin;

    /// The initial smoothing coefficients
    T m_coeff_init_lin;
    T m_coeff_init_NL;

    /// Do we print some messages
    bool m_print_messages;

    T m_ratio;

    /// Another function that can be called once the coefficient
    /// reduction iteration have converged
    gsFittingImproveAfterCV<d, GenBasis, T>* m_other_call;
};


}// namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

/** @file gsFittingParam.h

    @brief Contains the class containing all the parameters used
    by the fitting algorithm

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <gsIO/gsCmdLine.h>


namespace gismo
{


/**
   @brief Class containing the set of parameters used for the fitting
**/
template<class T = real_t>
class gsFittingParam
{
public:
    /// Initialization of all the parameters to their default value
    gsFittingParam();

    /// Copies the parameters of param
    gsFittingParam(gsFittingParam<T>& param);

    /// Initialization of the command line so that
    /// all the elements of this class can be set by the user
    void initOptEntry(gsCmdLine& cmd);

    /// Prints all the parameters
    void print(std::ostream& out);

    /// Checks that all the entry are conform
    void checkEntry();

    /// ATTRIBUTES
public:
    /// We stop the iterative process once the norm
    /// of the displacement is lower than this threshold
    T threshold_NL;

    /// A pointor to the identity, in case it is given
    /// by the user (template domain)
    gsFunctionSet<T>* p_identity;

    /// A pointor to the initialization, in case it is given
    /// by the user
    gsFunctionSet<T>* p_initialization;

    /// the prefix for all output files
    std::string output;

    /// In the case where the metric is considered for the minimization:
    /// consider the metric of the deformation if true
    bool deform_min;

    /// Identity used for:
    /// - minimization, when we consider the deformation
    /// - initialization of the iterative algorithm
    ///   in the nonlinear case
    /// If given (multipatch case most of the time),
    /// do we use the template domain as identity
    /// or do we construct another mapping?
    bool use_template;

    /// In the case where the template and target geometries
    /// are given:
    /// - if true, use directly these geometries in the minimization
    /// - if false, perform a sampling of the two objects and then
    ///   apply the least square method
    bool continuous_fitting;

    /// If true, the iterative refinement will be used
    bool use_refinement;

    /// the importance of the linear energy in the minimization
    T coeff_linear_global;

    /// the importance of the NL energy in the minimization
    T coeff_NL_global;

    /// coeff_gradient + coeff_hessian ~ 1
    /// (the importance of the linear energy
    ///  is given by coeff_linear_smoo)

    /// the importance of the gradient in the minimization
    T coeff_linear_gradient;

    /// the importance of the Hessian in the minimization
    T coeff_linear_hessian;

    /// coeff_NL_metric + coeff_NL_winslow ~ 1
    /// (the importance of the NL energy is given
    ///  by coeff_NL_smoo)

    /// the importance of the metric in the minimization
    T coeff_NL_metric;

    /// the importance of the Winslow energy in the minimization
    T coeff_NL_winslow;

    /// The domain must be refined till this threshold is reached
    /// (only used for the iterative refinement)
    T threshold_refinement;

    /// Proportion of refinement
    /// (only used for the iterative refinement)
    T proportion_refinement;

    /// We split elementary cells in each direction onto
    /// (extension_refinement + 1) parts
    /// (only used for the iterative refinement)
    index_t extension_refinement;

    /// Maximal number of refinement iterations
    index_t max_num_iter_refin;

    /// The maximal number of iterations for nonlinear energies
    index_t max_num_iter_NL;

    /// Export the points of the least squares method if true
    bool export_points;

    /// Export the solution found at each iteration if true
    bool export_iterations;

    /// Export the initialization if true
    bool export_initialization;

    /// Export the nodes after each refinement if true
    bool export_nodes;

    /// Should we print messages during the process
    bool print_messages;

    /// In case a basis is entered by the user, do we keep this
    /// basis unchanged
    bool keep_basis;

    /// Main input file. It can contain, depending on
    /// the problem considered:
    /// - the points of the template and of the target geometry
    /// - a spline representing the target geometry
    std::string input_file;

    /// In case the input_file contains the target geometry,
    /// the template geometry is given by the template_file
    std::string template_file;

    /// In case where we parametrize a multipatch domain,
    /// its topology is given in topo_file
    std::string topo_file;

    /// The degree of the splines used for the splitting
    index_t degree;

    /// The number of interior knots
    index_t interiorKnots;

    /// In case we first have to parametrize the boundary,
    /// the number of interiorKnots for this parametrization
    index_t interiorKnots_boundary;

    /// In the case of a nonlinear energy,
    /// if the norm of the first displacement is above this threshold,
    /// we recompute the initialization of the process
    T threshold_reinit;

    /// In the case where we recompute the initialization, the
    /// target points are simplified by interpolating the border of
    /// the current initialization with the target points.
    /// Coefficient of this interpolation
    T prop_reinit;

    /// In some cases, the fitting is called recursively.
    /// Index of the recursion.
    int ind;

    /// The number of slices used for to export 3D objects
    int nSlices;

    /// If true, we apply an algorithm that consists in projecting the solution
    /// after each computation. Then, the difference between the
    /// solution before and after the projection is added, up to a
    /// coefficient, to Lagrange multipliers
    bool use_lagrangeM;

    /// If true, we reduce the smoothing coefficient iteratively
    bool use_reducing_smoo;

    /// The coefficient used in the algorithm that consists in computing
    /// the Lagrange multipliers
    T coeff_lagrangeM;

    /// After the iterative method converge, we multiply the smoothing
    /// by this coefficient.
    T coeff_reducing_smoothing;

    /// The minimal smoothing coefficient.
    /// Once this threshold is reached, we do not reduce the
    /// smoothing parameter anymore.
    T inf_smoothing_parameter;

    /// The maximal number of smoothing parameter reduction autorized
    index_t max_num_reducing;

    /// The number of parameter lines used for the exportation
    /// of each iteration
    index_t num_lines;

    /// The coefficient associated with the Tikhonov regularization
    T coeff_tikhonov;

    /// If the projection algorithm is applied,
    /// frequency of its application. If one, it is applied
    /// every single iteration. If two, it is applied one over
    /// two iterations...
    index_t frequ_proj;

    /// The boundaries that are fixed. Each fixed boundary is
    /// given by
    /// - a triplet of int:
    /// the patch of the template that is fixed, the side of the template that is fixed,
    /// the patch of the target that is used
    /// - Followed by the beginning and the end of each variable (d-1 variables)
    std::vector<index_t> fixed_boundaries;

}; /// class gsFittingParam


}// namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFittingParam.hpp)
#endif

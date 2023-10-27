 /** @file gsFittingParam.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson

*/
#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsParameterLines.h>

#include <gsFitting/gsFittingBase.h>
#include <gsFitting/gsFittingQuadrature.hpp>
#include <gsFitting/gsFittingIdConstr.hpp>
#include <gsFitting/gsLeastSquares.hpp>


namespace gismo
{

template<class T> gsFittingParam<T>::
gsFittingParam(gsFittingParam<T>& param)
: output(param.output),
  input_file(param.input_file),
  template_file(param.template_file),
  topo_file(param.topo_file),
  fixed_boundaries(param.fixed_boundaries)
{
    p_identity = param.p_identity;
    p_initialization = param.p_initialization;

    continuous_fitting = param.continuous_fitting;
    deform_min = param.deform_min;
    use_template = param.use_template;

    use_refinement = param.use_refinement;

    coeff_linear_global = param.coeff_linear_global;
    coeff_NL_global = param.coeff_NL_global;

    coeff_linear_gradient = param.coeff_linear_gradient;
    coeff_linear_hessian = param.coeff_linear_hessian;

    coeff_NL_winslow = param.coeff_NL_winslow;
    coeff_NL_metric = param.coeff_NL_metric;

    threshold_NL = param.threshold_NL;
    max_num_iter_NL = param.max_num_iter_NL;

    threshold_refinement = param.threshold_refinement;
    proportion_refinement = param.proportion_refinement;
    extension_refinement = param.extension_refinement;
    max_num_iter_refin = param.max_num_iter_refin;

    export_nodes = param.export_nodes;
    export_initialization = param.export_initialization;
    export_iterations = param.export_iterations;
    export_points = param.export_points;
    print_messages = param.print_messages;

    keep_basis = param.keep_basis;

    degree = param.degree;
    interiorKnots = param.interiorKnots;
    interiorKnots_boundary = param.interiorKnots_boundary;

    prop_reinit = param.prop_reinit;
    threshold_reinit = param.threshold_reinit;

    ind = param.ind;

    use_lagrangeM = param.use_lagrangeM;
    coeff_lagrangeM = param.coeff_lagrangeM;

    nSlices = param.nSlices;

    num_lines = param.num_lines;

    coeff_reducing_smoothing = param.coeff_reducing_smoothing;
    inf_smoothing_parameter = param.inf_smoothing_parameter;
    max_num_reducing = param.max_num_reducing;
    coeff_tikhonov = param.coeff_tikhonov;

    use_reducing_smoo = param.use_reducing_smoo;
    frequ_proj = param.frequ_proj;
}

template<class T> gsFittingParam<T>::
gsFittingParam() : output("result"), input_file(""),
                   template_file(""), topo_file(""),
                   fixed_boundaries(0)
{
    p_identity = NULL;
    p_initialization = NULL;

    continuous_fitting = false;
    deform_min = true;
    use_template = false;

    use_refinement = false;

    coeff_linear_global = 1e-6;
    coeff_NL_global = 1e-6;

    coeff_linear_gradient = 0.1;
    coeff_linear_hessian = 0.9;

    coeff_NL_winslow = 1.;
    coeff_NL_metric = 0.;

    threshold_NL = 0.01;
    max_num_iter_NL = 10;

    threshold_refinement = 1e-6;
    proportion_refinement = 1.;
    extension_refinement = 1;
    max_num_iter_refin = 6;

    export_nodes = false;
    export_initialization = false;
    export_iterations = false;
    export_points = false;
    print_messages = false;

    keep_basis = true;

    degree = 2;
    interiorKnots = 10;
    interiorKnots_boundary = 30;

    prop_reinit = 0.5;
    threshold_reinit = 100000.;

    ind = 0;

    use_lagrangeM = false;
    coeff_lagrangeM = 0.1;

    nSlices = 5;

    num_lines = 30;

    coeff_reducing_smoothing = 0.4;
    inf_smoothing_parameter = 1e-6;
    max_num_reducing = 10;

    use_reducing_smoo = true;
    coeff_tikhonov = 0.;

    frequ_proj = 2;
}


template<class T>
void gsFittingParam<T>::print(std::ostream& out)
{
    checkEntry();
    out << "Input arguments: \n"
        << "output path:                     "
        << output << "\n"
        << "Main input:                      "
        << input_file << "\n"
        << "template file:                   "
        << template_file << "\n"
        << "Topology file: " << topo_file << "\n"
        << "continuous_fitting:              "
        << continuous_fitting << "\n"
        << "minimization on the deformation: "
        << deform_min << "\n"
        << "use template:                    "
        << use_template << "\n"
        << "coeff_linear_global:             "
        << coeff_linear_global << "\n"
        << "coeff_NL_global:                 "
        << coeff_NL_global << "\n"
        << "coeff_linear_gradient:           "
        << coeff_linear_gradient << "\n"
        << "coeff_linear_hessian:            "
        << coeff_linear_hessian << "\n"
        << "coeff_NL_winslow:                "
        << coeff_NL_winslow << "\n"
        << "coeff_NL_metric:                 "
        << coeff_NL_metric << "\n"
        << "threshold_NL:                    "
        << threshold_NL << "\n"
        << "max_num_iter_NL:                 "
        << max_num_iter_NL << "\n"
        << "export each iteration:           "
        << export_iterations << "\n"
        << "degree:                          "
        << degree << "\n"
        << "interior knots:                  "
        << interiorKnots << "\n"
        << "interiorKnots_boundary:          "
        << interiorKnots_boundary << "\n"
        /*  << "prop_reinit:                     "
        << prop_reinit << "\n"
        << "threshold_reinit:                "
        << threshold_reinit << "\n"   */
        << "keep_basis:                      "
        << keep_basis << "\n"
        << "use_lagrangeM:                   "
        << use_lagrangeM << "\n"
        << "coeff_lagrangeM:                 "
        << coeff_lagrangeM << "\n"
        << "max_num_reducing:                "
        << max_num_reducing << "\n"
        << "coeff_reducing_smoothing:        "
        << coeff_reducing_smoothing << "\n"
        << "inf_smoothing_parameter:         "
        << inf_smoothing_parameter << "\n";

    if(use_refinement)
    {
        out << "--------ITERATIVE REFINEMENT USED---------" << std::endl;
        out << "threshold_refinement:            "
            << threshold_refinement << "\n"
            << "proportion_refinement:           "
            << proportion_refinement << "\n"
            << "extension_refinement:            "
            << extension_refinement << "\n"
            << "max_num_iter_refin:              "
            << max_num_iter_refin << "\n";
    }
}

template<class T>
void gsFittingParam<T>::checkEntry()
{
    GISMO_ENSURE(coeff_linear_global >= 0, "coefficients must be positive");
    GISMO_ENSURE(coeff_linear_hessian >= 0, "coefficients must be positive");
    GISMO_ENSURE(coeff_linear_gradient >= 0, "coefficients must be positive");
    GISMO_ENSURE(coeff_NL_global >= 0, "coefficients must be positive");
    GISMO_ENSURE(coeff_NL_metric >= 0, "coefficients must be positive");
    GISMO_ENSURE(coeff_NL_winslow >= 0, "coefficients must be positive");

    GISMO_ENSURE(coeff_linear_global == 0
                 || coeff_linear_hessian > 0
                 || coeff_linear_gradient > 0,
                 "global linear coefficient positive with all the linear coefficients null");
    GISMO_ENSURE(coeff_NL_global == 0
                 || coeff_NL_metric > 0
                 || coeff_NL_winslow > 0,
                 "global nonlinear coefficient positive with all the nonlinear coefficients null");

    GISMO_ENSURE(threshold_refinement > 0, "thresholds must be strictly positive");
    GISMO_ENSURE(threshold_NL > 0, "thresholds must be strictly positive");

}


template<class T>
void gsFittingParam<T>::initOptEntry(gsCmdLine& cmd)
{
    cmd.addSwitch("R", "use_refinement", "Do we use the iterative refinement? (false by default)",
                  use_refinement);
    cmd.addSwitch("f", "continuous_fitting", "Do we use a continuous algorithm or do we sample the boundary before fitting? (false by default)",
                  continuous_fitting);
    cmd.addSwitch("D", "deformation", "Do we apply the minimization on the mapping or on the deformation, i.e., the difference between the mapping and the template mapping (identity)? (true by default)",
                  deform_min);
    /// only used in the case of multipatch fitting
    /// whenever the template domain is given
    cmd.addSwitch("t", "use_template", "In case of a multi-patch fitting, the template must be given. Do we use this template domain as reference mapping (for initialization and computing the deformation) or do we construct another one? (false by default)", use_template);

    cmd.addSwitch("E", "export_iterations", "Export each iteration if true (false by default)", export_iterations);
    cmd.addReal("l", "coeff_linear", "Importance of the linear energy (10^-6 by default)",
                coeff_linear_global);
    cmd.addReal("L", "coeff_NL", "Importance of the nonlinear energy. (10^-6 by default)",
                coeff_NL_global);
    cmd.addReal("G", "coeff_grad", "Importance of the gradient in the linear energy (between 0 and 1). (0.1 by default)", coeff_linear_gradient);
    cmd.addReal("H", "coeff_hess", "Importance of the hessian in the linear energy (between 0 and 1). (0.9 by default)",
                coeff_linear_hessian);
    cmd.addReal("W", "coeff_winslow", "Importance of the Winslow energy in the nonlinear energy (between 0 and 1). (1. by default)",
                coeff_NL_winslow);
    cmd.addReal("M", "coeff_metric", "Importance of the metric deformation in the nonlinear energy (between 0 and 1). (0. by default)",
                coeff_NL_metric);

    cmd.addString("o", "output", "Prefix for all output files ('result' by default)", output);
    cmd.addString("", "template_file", "Name of the file containing the template geometry. Only used if the input file is a geometry containing the target domain (default '')", template_file);
    cmd.addString("", "input_file", "Name of the file containing the input (default '')", input_file);
    cmd.addString("", "topo_file", "Name the file containing the topology of the basis used for constructing the parameterization. Only necessary if using multipatch (default '')", topo_file);


    cmd.addReal("T", "threshold_NL", "We stop the iterative process once the norm of the displacement is lower than this threshold. (0.01 by default)",
                threshold_NL);
    cmd.addInt("S", "max_num_iter_NL", "The maximal number of iterations for nonlinear energies (10 by default)",
                max_num_iter_NL);

    cmd.addReal("e", "threshold_refinement", "Error threshold. The domain must be refined till this threshold is reached (10^-6 by default)",
                threshold_refinement);
    cmd.addReal("P", "proportion_refinement", "Proportion of cells above the threshold that will be refined (1. by default)",
                proportion_refinement);
    cmd.addInt("i", "extension_refinement", "We split elementary cells in each direction onto (extension_refinement + 1) parts (1 by default)", extension_refinement);
    cmd.addInt("r", "max_num_iter_refin", "The maximal number of refinement iterations (6 by default)",
               max_num_iter_refin);

    cmd.addInt("d", "degree", "The degree of the spline (2 by default)", degree);
    cmd.addInt("k", "interKnots", "The number of interior knots (10 by default)",
               interiorKnots);
    cmd.addInt("K", "interKnots_boundary", "The number of interior knots used to fit the boundary in the case where a boundary must firstly be fitted (30 by default)",
               interiorKnots_boundary);

    cmd.addSwitch("B", "keep_basis", "In case a basis is given, do we keep this basis unchanged (true by default)",
               keep_basis);

/*    cmd.addReal("", "prop_reinit", "In the case where we recompute the initialization, the target points are simplified by interpolating the border of the current initialization with the target points. Coefficient of this interpolation (0.5 by default)",
               prop_reinit);
               cmd.addReal("", "threshold_reinit", "In the case of a nonlinear energy, if the norm of the first displacement is above this threshold, we recompute the initialization of the process (3. by default)",
        threshold_reinit);*/

    cmd.addSwitch("", "use_lagrangeM", "We apply an algorithm that consists in projecting the solution after each computation. Then, the difference between the solution before and after the projection is added, up to a coefficient, to Lagrange multipliers. (false by default)",
                  use_lagrangeM);

    cmd.addReal("", "coeff_lagrangeM", "The coefficient used the for algorithm that consists in projecting the solution and then adapting the Lagrange multipliers (0.1 by default)",
                coeff_lagrangeM);

    cmd.addSwitch("", "print_messages", "Print all the messages of the iterations if true. (false by default)",
                  print_messages);
    cmd.addSwitch("", "export_initialization", "Export the initialization of the algorithm. (false by default)",
                  export_initialization);
    cmd.addSwitch("", "export_points", "Export the points used for the least squares method. (false by default)",
                  export_points);
    cmd.addSwitch("", "export_nodes", "Export the nodes after each refinement. (false by default)",
                  export_nodes);

    cmd.addInt("", "frequency_projection", "If the projection algorithm is applied, frequency of its application. If one, it is applied every single iteration. If two, it is applied one over two iterations... (2 by default)",
               frequ_proj);

    cmd.addInt("", "max_num_reducing", "After the algorithm converge, we reduce the smoothing parameter to improve the accuracy of the smoothing. The maximal number of smoothing parameter reduction autorized (10 by default)",
               max_num_reducing);
    cmd.addReal("", "coeff_reducing_smoothing", "After the algorithm converge, we reduce multiply recursively the smoothing parameter by this parameter to improve the accuracy of the smoothing. (0.4 by default)",
               coeff_reducing_smoothing);
    cmd.addReal("", "inf_smoothing_parameter", "Once this threshold is reached, we don't reduce smoothing coefficients anymore. (1e-6 by default)",
               inf_smoothing_parameter);
    cmd.addInt("", "num_lines", "The number of coordinate curves used for the exportation of each iteration (30 by default)",
               num_lines);
    cmd.addReal("", "coeff_tikhonov", "The coefficient associated with the Tikhonov regularization (30 by default)",
                coeff_tikhonov);
    cmd.addSwitch("", "use_reducing_smoo", "If true, the smoothing coefficients are reduced iteratively. (true by default)",
                  use_reducing_smoo);
    cmd.addMultiInt("", "fixed_boundaries", "The boundaries that are fixed. Each fixed boundary is first associated with a triplet of int: the patch of the template that is fixed, the side of the template that is fixed, the patch of the target that is used. This is followed by the beginning and the end of each variable (d-1 variables)",
                    fixed_boundaries);

}

} // namespace gismo

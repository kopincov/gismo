/** @file tutorial_opt.cpp

    @brief Constructs a Spline from input files (see the instructions below).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#include <gismo.h>
#include <gsFitting/gsFittingParam.hpp>
#include <gsFitting/gsFittingUtilsIO.h>
#include <gsFitting/gsFittingUtilsIO.hpp>

using namespace gismo;

int main(int argc, char* argv[])
{
    gsFittingParam<> fitting_param;
    gsCmdLine cmd("Constructs a spline (single patch or multi-patch) that fits some data in 2D or 3D. This data can be:\n - a boundary or a trimming curve (or surface), in which case the input file (given by the argument --input_file) should be a geometry (single or multipatch) called target geometry. This geometry can be multipatch even in the case where a single patch spline is researched. This can occur for example when multiple boundaries (or trimmed curves,...) must be fitted. On the other hand, if the expected spline is multipatch, then the patches must coincide (to be improved). The template geometry, i.e., the preimage of the target geometry, should be given with the argument --template_file.\n - a set of points to be fitted and there preimage in the parameter domain, both defined in one file given by the argument --input_file.\n - a mapping between two different solids that must be found. The input file (given by the argument --input_file) should contain the target solid (the image of the mapping) and the template file (given by the argument --template_file) should contain the template solid (the domain of the mapping).\n Note that, to consider multipatch fitting, the parameter domain (and its topology) should be given using --topo_file. This input can either be a multipatch geometry or a multi-basis.\n The fitting is obtained by minimizing some energy (obtained by the least squares method). To this energy is added a smoothing that permits both to ensure the convergence of the minimization and to improve the regularity of the result. This smoothing energy is given by a linear combinaison of:\n - linear energies with coefficent given by --coeff_linear_global: norm of the gradient (proportion in [0,1] given by --coeff_linear_gradient) and norm of the second derivatives (proportion in [0,1] given by --coeff_linear_hessian)\n - nonlinear energies with coefficent given by --coeff_NL_global: norm of the metric (proportion in [0,1] given by --coeff_NL_metric) and Winslow energy (proportion in [0,1] given by --coeff_NL_metric).\n A hierarchical refinement can also be applied (by the argument --use_refinement). Finally note that to fit a heavily deformed shape, a high smoothing coefficient should be chosen. This coefficient is then reduced once the convergence is reached. This permits to initialize the algorithm.\n REMARKS: things that should be modified:\n - Lagrange multipliers work fine in the case of a linear energy (though the convergence is quite slow). In the case of nonlinear energies, multiple iterations of multipliers adapting should be performed before displacing the position. Current implementation that displace the current position after one adaptation of Lagrange multipliers does not converge\n - Case of multipatch fitting with continuous fitting: must be checked.\n - Automatic increasing of the smoothing coefficient in case the displacement is too important, decreasing faster the smoothing coefficient in case the displacement is too small.\n - Make the parameters of the algorithm independent of the scaling.");

    /// Sets the default values
    fitting_param.initOptEntry(cmd);
    /// read the input data
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if(fitting_param.input_file.empty())
    {
        gsInfo << "Please provide a file name or call -h for help.\n";
        return 0;
    }
    
    fitting_param.print(gsInfo);

    /// computes the fitting
    readInput<real_t>(fitting_param);
    return 0;
}



/** @file gsOptADtest

    @brief Test of optimization with AD in reverse mode with use of CoDiPack and IpOpt
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Jaeschke
*/

#include <iostream>

#undef GISMO_BUILD_LIB
#include <gsCore/gsConfig.h>
#include <gsCore/gsDevConfig.h>

#include <gismo.h>

#if defined(GISMO_WITH_ADDSL)
#include <gsAdDSL/gsAdDSL.h>
#endif

// #include <gsModeling/gsSpringPatch.h>
// #include <gsModeling/gsCoonsPatch.h>
//#include <gsOptimizer/gsOptParamJacobian.h>
//#include <gsOptimizer/gsOptParamWinslow.h>
//#include <gsOptimizer/gsOptParamLiao.h>
//#include <gsOptimizer/gsOptParamAreaOrth.h>
#include <gsOptimizer/gsOptADArea.h>

using namespace gismo;

int main(int argc, char* argv[])
    {

	// Input options
	real_t penalty = 5.0;
	std::string path = "boundaries/square_bd.xml";
	bool info_dom = false;
	bool info_diff = false;
	bool plot = true;
	int numElevate  = 0;
    int numHref     = 1;
    int basisDegree = 0;
	gsCmdLine cmd("Testing an optimization problem.");
	cmd.addString("g","geometry",
		"Path to the initial geometry",
		path);
	cmd.addReal("o","penalty","Value of penalty parameter",penalty);
	cmd.addInt("r","hRefine",
		"Number of dyadic h-refinement (bisection) steps to perform before solving",
		numHref);
	cmd.addInt("p","degree",
		"Degree of the basis functions to use for solving (will elevate or reduce the input)",
		basisDegree);
	cmd.addInt("e","degreeElevation",
		"Number of degree elevation steps to perform on the Geometry's basis before solving",
		numElevate);
	cmd.addSwitch("plot", "Plot result in ParaView format", plot);
	cmd.addSwitch("info_dom", "Display additional information about the domain reconstruction", info_dom);
	cmd.addSwitch("info_diff", "Display additional information about the differentiation", info_diff);

	try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

	gsOptADArea a (plot,path,penalty,numElevate,numHref,basisDegree, info_dom, info_diff);

  gsInfo << "---------------------------------------------";
  a.solve();

  // Print some details in the output

    // Print final design info
    gsInfo << "Number of iterations : " << a.iterations() <<"\n";
    gsInfo << "Final objective value: " << a.objective() <<"\n";

    return 0;
}

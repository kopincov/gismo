/** @file utils_plot_points

    @brief Plots points from a file. 

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#include <gismo.h>


#include <iostream>
#include <string>
#include "gsFittingUtils.h"

using namespace gismo;



int main(int argc, char* argv[])
{
    std::string input("");
    std::string output("out");

    gsCmdLine cmd("Plotting points.");
    cmd.addPlainString("input", "Input file", input);
    cmd.addString("o", "output", "Output file", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << " \n\nInput arguments: \n\n"
	   << "input:  " << input << "\n\n"
	   << "output: " << output << "\n\n"
	   << "--------------------------------------------------\n" << "\n";
    
    if (input == "")
    {
	gsInfo << "No input\n\n";
    }

    gsFileData<> fd(input);
    const int num_patches = number_of_patches(input);
    
    for (int i = 0; i != num_patches; i++)
    {
	gsMatrix<> parameters = *fd.getId< gsMatrix<> >(2 * i);
	gsMatrix<> points = *fd.getId< gsMatrix<> >(2 * i + 1);

	std::string out = output + "_params_" + util::to_string(i);
	gsInfo << "Writing: " << out << "\n";
	gsWriteParaviewPoints(parameters, out);
    
	out = output + "_points_" + util::to_string(i);
	gsInfo << "Writing: " << out << "\n";
	gsWriteParaviewPoints(points, out);
    }

    return 0;
}



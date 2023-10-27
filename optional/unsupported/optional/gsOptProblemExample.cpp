/** @file gsOptProblemExample

    @brief Toy example for optimizer

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

#include <gsOptimizer/gsOptProblemExample.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    gsOptProblemExample<real_t> optimizer;
            
    // Run optimizer
    optimizer.solve();

    // Print some details in the output
    optimizer.print(gsInfo);

    // Print final design info
    gsInfo << "Number of iterations : " << optimizer.iterations() <<"\n";
    gsInfo << "Final objective value: " << optimizer.objective() <<"\n";

    return 0;
}

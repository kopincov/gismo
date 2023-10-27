/** @file gsCompareResults.cpp

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char** argv)
{
    /************** Define command line options *************/

    std::string fn1;
    std::string fn2;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addString("a",  "fn1",          "First xml file", fn1);
    cmd.addString("b",  "fn2",          "Second xml file", fn2);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (fn1.empty() || fn2.empty()) { gsInfo << "Use --fn1 and --fn2 to specify files"; return 0; }

    memory::unique_ptr< gsMatrix<> > m1 = gsReadFile<>(fn1);
    if (!m1) return EXIT_SUCCESS;
    memory::unique_ptr< gsMatrix<> > m2 = gsReadFile<>(fn2);
    if (!m2) return EXIT_SUCCESS;

    gsInfo << "Matrix 1:        " << m1->rows() << " x " << m1->cols() << "\n";
    gsInfo << "Matrix 2:        " << m2->rows() << " x " << m2->cols() << "\n";

    if (m1->rows() != m2->rows() || m1->cols() != m2->cols()) return EXIT_SUCCESS;

    gsInfo << "Difference:      " << (*m1-*m2).norm() << "\n";
    gsInfo << "Rel. difference: " << (*m1-*m2).norm()/math::min(m1->norm(),m2->norm()) << "\n";

    return EXIT_SUCCESS;
}

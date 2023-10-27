/** @file gsScrewMachineCFX5.cpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moeller
*/

#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;

#ifdef GISMO_WITH_CFX5

#include "gsScrewMachineCFX5.hpp"

int main(int argc, char* argv[])
{  
  gsCFX5ScrewMachine cfx5(argc, argv);
  return 0;
}

#else

int main(int argc, char *argv[])
{
  gsInfo << "gsScrewMachineCFX5 must be compiled with GISMO_WITH_CFX5=ON\n";
  return 1;
}

#endif

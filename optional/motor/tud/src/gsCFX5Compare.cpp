/** @file gsCFX5Compare.cpp

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

#include <gsCFX5/gsCFX5.h>

int main(int argc, char* argv[])
{
  std::string fn = "", fn_compare = "";
  std::string boundary = "";
  std::string domain   = "";
  std::string region   = "";
  std::string variable = "";
  std::string volume   = "";

  gsCmdLine cmd("gsCFX5Compare: Comparison of CFX5 mesh and result files.");
  cmd.addPlainString("filename", "Name of the CFX5 file.", fn);
  cmd.addString("c", "compare",  "Name of the comparison file.", fn_compare);
  cmd.addString("V", "variable", "Show single variable.",  variable);
  cmd.addString("b", "boundary", "Show single boundary.",  boundary);
  cmd.addString("d", "domain",   "Show single domain.",    domain);
  cmd.addString("r", "region",   "Show single region.",    region);
  cmd.addString("v", "volume",   "Show single volume.",    volume);

  // Get command line arguments
  try { cmd.getValues(argc,argv); } catch (int rv) { exit(0); }

  gsCFX5 cfx5, cfx5_compare;

  GISMO_ENSURE(!fn.empty(), "No input file given.");
  cfx5.importFromCFX5(fn);

  GISMO_ENSURE(!fn_compare.empty(), "No comparison file given.");
  cfx5_compare.importFromCFX5(fn_compare);

    if (domain != "")
    {
      // Get pointer to selected domain
      gsCFX5Domain *cfx5Domain, *cfx5Domain_compare;
      try
        {
          cfx5Domain = cfx5.searchDomain(std::stoi(domain));
          cfx5Domain_compare = cfx5_compare.searchDomain(std::stoi(domain));
        }
      catch (...)
        {
          cfx5Domain = cfx5.searchDomain(domain);
          cfx5Domain_compare = cfx5_compare.searchDomain(domain);
        }

      GISMO_ENSURE(cfx5Domain && cfx5Domain_compare, "Invalid domain: "+domain);

      if (region != "")
        {
          gsCFX5Region *cfx5Region, *cfx5Region_compare;
          try
            {
              cfx5Region = cfx5Domain->searchRegion(std::stoi(region));
              cfx5Region_compare = cfx5Domain_compare->searchRegion(std::stoi(region));
            }
          catch (...)
            {
              cfx5Region = cfx5Domain->searchRegion(region);
              cfx5Region_compare = cfx5Domain_compare->searchRegion(region);
            }

          GISMO_ENSURE(cfx5Region && cfx5Region_compare, "Invalid region: "+region);

          gsInfo << cfx5Region->diff(*cfx5Region_compare);
        }
      if (volume != "")
        {
          gsCFX5Volume *cfx5Volume, *cfx5Volume_compare;
          try
            {
              cfx5Volume = cfx5Domain->searchVolume(std::stoi(volume));
              cfx5Volume_compare = cfx5Domain_compare->searchVolume(std::stoi(volume));
            }
          catch (...)
            {
              cfx5Volume = cfx5Domain->searchVolume(volume);
              cfx5Volume_compare = cfx5Domain_compare->searchVolume(volume);
            }

          GISMO_ENSURE(cfx5Volume && cfx5Volume_compare, "Invalid volume: "+volume);

          gsInfo << cfx5Volume->diff(*cfx5Volume_compare);
        }
      if (boundary != "")
        {
          gsCFX5Boundary *cfx5Boundary, *cfx5Boundary_compare;
          try
            {
              cfx5Boundary = cfx5Domain->searchBoundary(std::stoi(boundary));
              cfx5Boundary_compare = cfx5Domain_compare->searchBoundary(std::stoi(boundary));
            }
          catch (...)
            {
              cfx5Boundary = cfx5Domain->searchBoundary(boundary);
              cfx5Boundary_compare = cfx5Domain_compare->searchBoundary(boundary);
            }

          GISMO_ENSURE(cfx5Boundary && cfx5Boundary_compare, "Invalid boundary: "+boundary);

          gsInfo << cfx5Boundary->diff(*cfx5Boundary_compare);
        }
      if (variable != "")
        {
          gsCFX5Variable *cfx5Variable, *cfx5Variable_compare;
          try
            {
              cfx5Variable = cfx5Domain->searchVariable(stoi(variable));
              cfx5Variable_compare = cfx5Domain_compare->searchVariable(stoi(variable));
            }
          catch (...)
            {
              cfx5Variable = cfx5Domain->searchVariable(variable);
              cfx5Variable_compare = cfx5Domain_compare->searchVariable(variable);
            }

          GISMO_ENSURE(cfx5Variable && cfx5Variable_compare, "Invalid variable: "+variable);

          gsInfo << cfx5Variable->diff(*cfx5Variable_compare);
        }
      if (region+volume+boundary+variable == "")
        {
          gsInfo << cfx5Domain->diff(*cfx5Domain_compare);
        }

    }
  else
    {
      // Complete comparison
      gsInfo << cfx5.diff(cfx5_compare);
    }

  return 0;
}

#else

int main(int argc, char *argv[])
{
  gsInfo << "gsCFX5Compare must be compiled with GISMO_WITH_CFX5=ON\n";
  return 1;
}

#endif

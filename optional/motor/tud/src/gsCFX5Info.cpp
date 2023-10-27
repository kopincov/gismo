/** @file gsCFX5Info.cpp

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
  bool verbose = false;
  std::string boundary = "";
  std::string domain   = "";
  std::string fn = "";
  std::string region   = "";
  std::string variable = "";
  std::string volume   = "";

  gsCmdLine cmd("gsCFX5Info: Information about CFX5 mesh and result file.");
  cmd.addPlainString("filename", "Name of the CFX5 file.", fn);
  cmd.addString("V", "variable", "Show single variable.",  variable);
  cmd.addString("b", "boundary", "Show single boundary.",  boundary);
  cmd.addString("d", "domain",   "Show single domain.",    domain);
  cmd.addString("r", "region",   "Show single region.",    region);
  cmd.addString("v", "volume",   "Show single volume.",    volume);
  cmd.addSwitch(     "verbose",  "Verbose output.",        verbose);

  // Get command line arguments
  try { cmd.getValues(argc,argv); } catch (int rv) { exit(0); }

  gsCFX5 cfx5;

  GISMO_ENSURE(!fn.empty(), "No input file given.");
  cfx5.importFromCFX5(fn);

  if (domain != "")
    {
      // Get pointer to selected domain
      gsCFX5Domain *cfx5Domain;
      try
        {
          cfx5Domain = cfx5.searchDomain(std::stoi(domain));
        }
      catch (...)
        {
          cfx5Domain = cfx5.searchDomain(domain);
        }

      GISMO_ENSURE(cfx5Domain, "Invalid domain: "+domain);

      if (region != "")
        {
          gsCFX5Region *cfx5Region;
          try
            {
              cfx5Region = cfx5Domain->searchRegion(std::stoi(region));
            }
          catch (...)
            {
              cfx5Region = cfx5Domain->searchRegion(region);
            }

          GISMO_ENSURE(cfx5Region, "Invalid region: "+region);

          if (verbose)
            gsInfo << cfx5Region->detail();
          else
            cfx5Region->print(gsInfo);
          gsInfo << std::endl;
        }
      if (volume != "")
        {
          gsCFX5Volume *cfx5Volume;
          try
            {
              cfx5Volume = cfx5Domain->searchVolume(std::stoi(volume));
            }
          catch (...)
            {
              cfx5Volume = cfx5Domain->searchVolume(volume);
            }

          GISMO_ENSURE(cfx5Volume, "Invalid volume: "+volume);

          if (verbose)
            gsInfo << cfx5Volume->detail();
          else
            cfx5Volume->print(gsInfo);
          gsInfo << std::endl;
        }
      if (boundary != "")
        {
          gsCFX5Boundary *cfx5Boundary;
          try
            {
              cfx5Boundary = cfx5Domain->searchBoundary(std::stoi(boundary));
            }
          catch (...)
            {
              cfx5Boundary = cfx5Domain->searchBoundary(boundary);
            }

          GISMO_ENSURE(cfx5Boundary, "Invalid boundary: "+boundary);

          if (verbose)
            gsInfo << cfx5Boundary->detail();
          else
            cfx5Boundary->print(gsInfo);
          gsInfo << std::endl;
        }
      if (variable != "")
        {
          gsCFX5Variable *cfx5Variable;
          try
            {
              cfx5Variable = cfx5Domain->searchVariable(stoi(variable));
            }
          catch (...)
            {
              cfx5Variable = cfx5Domain->searchVariable(variable);
            }

          GISMO_ENSURE(cfx5Variable, "Invalid variable: "+variable);

          if (verbose)
            gsInfo << cfx5Variable->detail();
          else
            cfx5Variable->print(gsInfo);
          gsInfo << std::endl;
        }
      if (region+volume+boundary+variable == "")
        {
          if (verbose)
            gsInfo << cfx5Domain->detail();
          else
            cfx5Domain->print(gsInfo);
          gsInfo << std::endl;
        }

    }
  else
    {
      // Show complete information
      if (verbose)
        gsInfo << cfx5.detail();
      else
        cfx5.print(gsInfo);
      gsInfo << std::endl;
    }

  return 0;
}

#else

int main(int argc, char *argv[])
{
  gsInfo << "gsCFX5Info must be compiled with GISMO_WITH_CFX5=ON\n";
  return 1;
}

#endif

/** @file gsRenderer.cpp

    @brief Vizualize G+Smo objects from XML input

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>
#include <gismo_dev.h>

using namespace gismo;

#ifdef GISMO_WITH_RENDERER

#include <gsRender.h>

int main(int argc, char *argv[])
{
  std::string fn("");

  //! [Parse Command line]
  gsCmdLine cmd("Hi, give me an XML file and I will try to visualize it!");

  cmd.addPlainString("filename", "File containing data to visualize (.xml)", fn);
  
  cmd.getValues(argc,argv);
  //! [Parse Command line]

  if ( fn.empty() )
    {
      gsInfo<< cmd.getMessage();
      gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
      return 0;
    }
  
  gsRender renderer;
  gsFileData<>  filedata(fn);

  if ( filedata.has< gsMultiPatch<> >() )
    {
      std::vector< gsMultiPatch<>::uPtr > mp = filedata.getAll< gsMultiPatch<> >();
      for (auto it = mp.begin(); it != mp.end(); it++)
        renderer.addObject(**it, fn);
    }
  else if ( filedata.has< gsGeometry<> >() )
    {
      std::vector< gsGeometry<>::uPtr > geo = filedata.getAll< gsGeometry<> >();
      gsMultiPatch<> mp;
  
      for (auto it = geo.begin(); it != geo.end(); ++it)
        {
          mp.addPatch(**it);
        }
      geo.clear();
      renderer.addObject(mp, fn);
    }
  else
    {
      gsInfo<< "Did not find anything to plot in " << fn << ", quitting." << "\n";
    }
  filedata.clear();
  
  renderer.readSettings((std::string)getenv("PROJECT_PATH")+"/gsRenderSettings.ini");
  renderer.startGUI();
  
  return 0;
}

#else

int main(int argc, char *argv[])
{
  gsInfo << "gsRenderer must be compiled with GISMO_WITH_RENDERER=ON\n";
  return 1;
}

#endif

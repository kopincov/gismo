/** @file gsScrewMachineCFX5.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moeller
*/

#pragma once

#include <gismo.h>
#include <gsCFX5/gsCFX5.h>

namespace gismo {

/**
   \brief Screw machine CFX5/G+Smo import/export tool

   The screw machine multi-patch geometry consists of five patches
   with the following topology

   +-----+-----+-----+
   |     |     |     |
   |  1  |  2  |  3  |
   |     |     |     |
   +-----+-----+-----+
   |     |     |     |
   |  0  |     |  4  |
   |     |     |     |
   +-----+     +-----+

 */
class gsCFX5ScrewMachine : public gsCFX5
{
public:

  /// \brief Constructor
  gsCFX5ScrewMachine(int argc, char* argv[])
    : gsCFX5::gsCFX5()
  {
    /// Names of input/output files
    std::string fninput  = "";
    std::string fnoutput = "";
    std::string fnlog    = "";
    std::string fnparam  = "";
    std::string fnconfig = "";

    /// Operation modes
    bool xml_icem  = false;
    bool xml_cfx5  = false;
    bool xml_xml   = false;
    bool res_xml   = false;
    
    // Initialize default command line arguments
    gsCmdLine cmd = getCmdLine("Screw machine CFX5/G+Smo import/export tool");

    cmd.addSwitch("xml_icem",
                  "Convert G+Smo XML file to Ansys ICEM mesh file.",
                  xml_icem);

    cmd.addSwitch("xml_cfx5",
                  "Convert G+Smo XML file to Ansys CFX5 mesh file.",
                  xml_cfx5);

    cmd.addSwitch("xml_xml",
                  "Convert G+Smo XML file to G+Smo XML file.",
                  xml_xml);

    cmd.addSwitch("res_xml",
                  "Convert Ansys CFX5 result file to G+Smo XML file.",
                  res_xml);

    cmd.addString("i", "input",
                  "Name of the input file.",
                  fninput);

    cmd.addString("o", "output",
                  "Name of the output file.",
                  fnoutput);

    cmd.addString("l", "log",
                  "Name of the log file.",
                  fnlog);

    cmd.addString("p", "parameterization",
                  "Name of the parameterization file.",
                  fnparam);

    cmd.addPlainString("config",
                       "File containing configuration.", fnconfig);

    // Get command line arguments
    try { cmd.getValues(argc,argv); } catch (int rv) { exit(0); }

    // Read configuration from file (using command-line overrides)
    getOptionList(fnconfig);

    // Convert G+Smo XML file to Ansys ICEM mesh file
    if (xml_icem)
      {
        timer.restart();
        importFromXML(fninput);
        gsInfo << "Import from XML in time: " << timer.stop() << std::endl;

        timer.restart();
        exportToICEM(fnoutput);
        gsInfo << "Export to ICEM in time: " << timer.stop() << std::endl;
      }
    else if (xml_cfx5)
      {
        timer.restart();
        importFromXML(fninput);
        gsInfo << "Import from XML in time: " << timer.stop() << std::endl;

        timer.restart();
        exportToCFX5(fnoutput);
        gsInfo << "Export to CFX5 in time: " << timer.stop() << std::endl;
      }
    else if (xml_xml)
      {
        gsFileData<> fd;
        timer.restart();
        fd.add(importFromXML(fninput));
        gsInfo << "Import from XML in time: " << timer.stop() << std::endl;

        timer.restart();
        fd.dump(fnoutput);
        gsInfo << "Export to XML in time: " << timer.stop() << std::endl;
      }
    else if (res_xml)
      {
        timer.restart();
        importFromCFX5(fninput);
        gsInfo << "Import from CFX5 in time: " << timer.stop() << std::endl;

        timer.restart();
        exportToXML(fnoutput, fnparam, {"TEMP_FL1", "TTOT_FL1"});
        gsInfo << "Export to XML in time: " << timer.stop() << std::endl;
      }
    else
      GISMO_ERROR("Unsupported operation.");

    // Write configuration to log
    if (!fnlog.empty())
      {
        timer.restart();
        gsFileData<> fd;
        fd.add(setOptionList());
        fd.dump(fnlog);
        gsInfo << "Writing log file in time: " << timer.stop() << std::endl;
      }
  }

  /// \brief Imports G+Smo multi-patch from XML file and stores it in
  /// the internal CFX5 mesh import/export API data structures
  gsMultiPatch<> importFromXML(const std::string& filename)
  {
    // Read G+Smo XML multi-patch parameterization
    gsFileData<> fileData(filename);
    gsMultiPatch<> geo, grid;

    if (fileData.has< gsMultiPatch<> >())
      {
        fileData.getFirst< gsMultiPatch<> >(geo);
        geo.degreeElevate(numElevate);
        for (int r=0; r<numRefine-1; ++r)
          geo.uniformRefine();
      }
    else
      throw std::runtime_error("Input file doesn't have a multipatch geometry inside.");

    if (geo.parDim() == 2)
      {
        gsCFX5Domain* domain = createDomain();
        domain->setName("Fluid");
        grid = genGridFromGeo2d(geo);
        genCFX5NodesFromMultiPatch2d(grid, domain);
        genCFX5Elements(domain);
      }
    else if (geo.parDim() == 3)
      {
        gsCFX5Domain* domain = createDomain();
        domain->setName("Fluid");
        grid = genGridFromGeo3d(geo);
        genCFX5NodesFromMultiPatch3d(grid, domain);
        genCFX5Elements(domain);
      }
    else
      GISMO_ERROR("Unsupported parametric dimension!");

    return grid;
  }

  /// \brief Exports CFX5 mesh and results that are stored internally
  /// to G+Smo multi-patch XML file
  std::vector<gsMultiPatch<> > exportToXML(const std::string& filename="",
                                           const std::string& parfilename="",
                                           const std::list<std::string>& variablenames={},
                                           const std::string& domainname="Default Domain")
  {
    // Convert CFX5 mesh and selected solution variables to G+Smo
    // multi-patch parameterization using degree one basis functions
    std::vector<gsMultiPatch<> > mp = convertToXML(variablenames, domainname);

    // Project CFX5 solution onto extra G+Smo parameterization
    if (!parfilename.empty())
      {
        std::vector<gsMultiPatch<> > mp_par;
        gsMultiPatch<> par;
        gsFileData<> fd(parfilename);
        fd.getFirst(par);
        par.degreeElevate(numElevate);
        for (int r=0; r<numRefine-1; ++r)
          par.uniformRefine();

        for (auto it = mp.begin(); it != mp.end(); it++)
          {
            gsMultiPatch<> mp_interp;
            auto par_it = par.begin();
            for (auto mp_it = it->begin(); mp_it != it->end(); mp_it++, par_it++)
              {
                if ( ((*par_it)->parDim() == (*mp_it)->parDim()) &&
                     ((*par_it)->geoDim() == (*mp_it)->geoDim()) )
                  {
                    // Convert between compatible CFX5 and G+Smo data
                    gsGeometry<>::uPtr interp =
                      (*par_it)->basis().interpolateAtAnchors(
                      (*mp_it)->eval((*par_it)->basis().anchors()));
                    mp_interp.addPatch(*interp);
                  }
                else if ( ((*par_it)->parDim() == 2) &&
                          ((*par_it)->geoDim() == 2) &&
                          ((*mp_it)->parDim()  == 3) &&
                          ((*mp_it)->geoDim()  == 3) )
                  {
                    // Convert 3d CFX5 mesh to 2d G+Smo geometry
                    gsGeometry<>::uPtr geo=(*mp_it)->boundary(boundary::front);
                    geo->embed(2);

                    gsGeometry<>::uPtr interp =
                      (*par_it)->basis().interpolateAtAnchors(
                        geo->eval((*par_it)->basis().anchors()));
                    mp_interp.addPatch(*interp);
                  }
                else if ( ((*par_it)->parDim() == 2) &&
                          ((*par_it)->geoDim() == 2) &&
                          ((*mp_it)->parDim()  == 3) &&
                          ((*mp_it)->geoDim()  == 1) )
                  {
                    // Convert CFX5 scalar data in 3d to G+Smo data in 2d
                    gsGeometry<>::uPtr geo=(*mp_it)->boundary(boundary::front);
                    geo->embed(2);

                    gsGeometry<>::uPtr interp =
                      (*par_it)->basis().interpolateAtAnchors(
                        geo->eval((*par_it)->basis().anchors()));
                    interp->embed(1);
                    mp_interp.addPatch(*interp);
                  }
                else
                  GISMO_ERROR("Incompatible combination of CFX5/G+Smo.");
              }
            mp_par.push_back(mp_interp);
          }
        mp_par.front().computeTopology();
        mp_par.front().addAutoBoundaries();
        mp.swap(mp_par);
      }

    // Export G+Smo parameterization to XML file
    if (!filename.empty())
      {
        gsFileData<> fd;
        for (auto it = mp.cbegin(); it != mp.cend(); it++)
          fd.add(*it);
        fd.dump(filename);
      }
    
    return mp;
  }

  /// \brief Converts CFX5 mesh and (if requested) results that are
  /// stored internally to G+Smo multi-patch objects
  std::vector<gsMultiPatch<> > convertToXML(const std::list<std::string>& variablenames={},
                                            const std::string& domainname="Default Domain")
  {
    gsCFX5Domain* domain = searchDomain(domainname);
    if (domain == nullptr)
      GISMO_ERROR("Invalid domain name: " + domainname + ".");

    // Generate uniform knot vectors
    gsKnotVector<> KVZ(0, 1, nz,      2, 1, 1);
    gsKnotVector<> KVL(0, 1, nleft,   2, 1, 1);
    gsKnotVector<> KVM(0, 1, nmiddle, 2, 1, 1);
    gsKnotVector<> KVR(0, 1, nright,  2, 1, 1);
    gsKnotVector<> KVT(0, 1, ntop,    2, 1, 1);
    gsKnotVector<> KVB(0, 1, nbottom, 2, 1, 1);

    // Define type 3D-tensor-B-spline pointer
    typedef typename gsTensorBSpline<3>::uPtr TensorBSpline3Ptr;

    std::vector<gsMultiPatch<> > mp;
    std::vector<gsCFX5Variable*> variable;
    std::vector<gsMatrix<> > V;
    gsCFX5Volume* volume;
    gsMultiPatch<> grid,mpsol;
    gsMatrix<> C;
    int idx,i,j,k,elem;

    // Geometry
    mp.push_back(gsMultiPatch<>());

    // Additional variables
    for (auto it = variablenames.cbegin(); it != variablenames.cend(); it++)
      {
        mp.push_back(gsMultiPatch<>());
        variable.push_back(domain->searchVariable(*it, true));
        V.push_back(gsMatrix<>());
      }

    // Patch #0 (left,bottom):
    volume = domain->searchVolume(0);

    if ((2+nleft) * (2+nbottom) * (2+nz) != volume->getNodeNumber())
      GISMO_ERROR("Number of grid points does not match configuration.");
    C.resize((2+nleft) * (2+nbottom) * (2+nz), 3);
    for (std::size_t v=0; v<variable.size(); v++)
      V[v].resize((2+nleft) * (2+nbottom) * (2+nz), 1);

    elem=0;
    for (auto elemit = volume->element_cbegin(); elemit != volume->element_cend(); elemit++, elem++)
      {
        i = elem%(1+nleft);
        j = elem/(1+nleft);
        k = elem/((1+nleft)*(1+nbottom));

        idx = k   *(2+nleft)*(2+nbottom)+ j   *(2+nleft)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(0)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(0)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(0)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(0)-1);

        idx = k   *(2+nleft)*(2+nbottom)+ j   *(2+nleft)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(1)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(1)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(1)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(1)-1);

        idx = k   *(2+nleft)*(2+nbottom)+(j+1)*(2+nleft)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(2)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(2)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(2)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(2)-1);

        idx = k   *(2+nleft)*(2+nbottom)+(j+1)*(2+nleft)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(3)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(3)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(3)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(3)-1);

        idx = (k+1)*(2+nleft)*(2+nbottom)+ j   *(2+nleft)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(4)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(4)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(4)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(4)-1);

        idx = (k+1)*(2+nleft)*(2+nbottom)+ j   *(2+nleft)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(5)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(5)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(5)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(5)-1);

        idx = (k+1)*(2+nleft)*(2+nbottom)+(j+1)*(2+nleft)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(6)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(6)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(6)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(6)-1);

        idx = (k+1)*(2+nleft)*(2+nbottom)+(j+1)*(2+nleft)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(7)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(7)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(7)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(7)-1);
      }
    mp.front().addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVL, KVB, KVZ, give(C))));
    for (std::size_t v=0; v<variable.size(); v++)
      mp[v+1].addPatch( mp.front()[0].basis().makeGeometry( give(V[v]) ) );

    // Patch #1 (left,top):
    volume = domain->searchVolume(1);

    if ((2+nleft) * (2+ntop) * (2+nz) != volume->getNodeNumber())
      GISMO_ERROR("Number of grid points does not match configuration.");
    C.resize((2+nleft) * (2+ntop) * (2+nz), 3);
    for (std::size_t v=0; v<variable.size(); v++)
      V[v].resize((2+nleft) * (2+ntop) * (2+nz), 1);

    elem=0;
    for (auto elemit = volume->element_cbegin(); elemit != volume->element_cend(); elemit++, elem++)
      {
        i = elem%(1+nleft);
        j = elem/(1+nleft);
        k = elem/((1+nleft)*(1+ntop));

        idx = k   *(2+nleft)*(2+ntop)+ j   *(2+nleft)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(0)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(0)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(0)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(0)-1);

        idx = k   *(2+nleft)*(2+ntop)+ j   *(2+nleft)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(1)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(1)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(1)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(1)-1);

        idx = k   *(2+nleft)*(2+ntop)+(j+1)*(2+nleft)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(2)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(2)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(2)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(2)-1);

        idx = k   *(2+nleft)*(2+ntop)+(j+1)*(2+nleft)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(3)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(3)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(3)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(3)-1);

        idx = (k+1)*(2+nleft)*(2+ntop)+ j   *(2+nleft)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(4)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(4)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(4)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(4)-1);

        idx = (k+1)*(2+nleft)*(2+ntop)+ j   *(2+nleft)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(5)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(5)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(5)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(5)-1);

        idx = (k+1)*(2+nleft)*(2+ntop)+(j+1)*(2+nleft)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(6)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(6)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(6)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(6)-1);

        idx = (k+1)*(2+nleft)*(2+ntop)+(j+1)*(2+nleft)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(7)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(7)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(7)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(7)-1);
      }
    mp.front().addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVL, KVT, KVZ, give(C))));
    for (std::size_t v=0; v<variable.size(); v++)
      mp[v+1].addPatch( mp.front()[1].basis().makeGeometry( give(V[v]) ) );

    // Patch #2 (middle,top):
    volume = domain->searchVolume(2);

    if ((2+nmiddle) * (2+ntop) * (2+nz) != volume->getNodeNumber())
      GISMO_ERROR("Number of grid points does not match configuration.");
    C.resize((2+nmiddle) * (2+ntop) * (2+nz), 3);
    for (std::size_t v=0; v<variable.size(); v++)
      V[v].resize((2+nmiddle) * (2+ntop) * (2+nz), 1);

    elem=0;
    for (auto elemit = volume->element_cbegin(); elemit != volume->element_cend(); elemit++, elem++)
      {
        i = elem%(1+nmiddle);
        j = elem/(1+nmiddle);
        k = elem/((1+nmiddle)*(1+ntop));

        idx = k   *(2+nmiddle)*(2+ntop)+ j   *(2+nmiddle)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(0)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(0)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(0)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(0)-1);

        idx = k   *(2+nmiddle)*(2+ntop)+ j   *(2+nmiddle)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(1)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(1)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(1)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(1)-1);

        idx = k   *(2+nmiddle)*(2+ntop)+(j+1)*(2+nmiddle)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(2)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(2)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(2)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(2)-1);

        idx = k   *(2+nmiddle)*(2+ntop)+(j+1)*(2+nmiddle)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(3)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(3)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(3)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(3)-1);

        idx = (k+1)*(2+nmiddle)*(2+ntop)+ j   *(2+nmiddle)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(4)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(4)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(4)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(4)-1);

        idx = (k+1)*(2+nmiddle)*(2+ntop)+ j   *(2+nmiddle)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(5)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(5)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(5)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(5)-1);

        idx = (k+1)*(2+nmiddle)*(2+ntop)+(j+1)*(2+nmiddle)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(6)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(6)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(6)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(6)-1);

        idx = (k+1)*(2+nmiddle)*(2+ntop)+(j+1)*(2+nmiddle)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(7)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(7)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(7)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(7)-1);
      }
    mp.front().addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVM, KVT, KVZ, give(C))));
    for (std::size_t v=0; v<variable.size(); v++)
      mp[v+1].addPatch( mp.front()[2].basis().makeGeometry( give(V[v]) ) );

    // Patch #3 (right,top):
    volume = domain->searchVolume(3);

    if ((2+nright) * (2+ntop) * (2+nz) != volume->getNodeNumber())
      GISMO_ERROR("Number of grid points does not match configuration.");
    C.resize((2+nright) * (2+ntop) * (2+nz), 3);
    for (std::size_t v=0; v<variable.size(); v++)
      V[v].resize((2+nright) * (2+ntop) * (2+nz), 1);

    elem=0;
    for (auto elemit = volume->element_cbegin(); elemit != volume->element_cend(); elemit++, elem++)
      {
        i = elem%(1+nright);
        j = elem/(1+nright);
        k = elem/((1+nright)*(1+ntop));

        idx = k   *(2+nright)*(2+ntop)+ j   *(2+nright)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(0)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(0)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(0)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(0)-1);

        idx = k   *(2+nright)*(2+ntop)+ j   *(2+nright)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(1)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(1)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(1)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(1)-1);

        idx = k   *(2+nright)*(2+ntop)+(j+1)*(2+nright)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(2)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(2)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(2)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(2)-1);

        idx = k   *(2+nright)*(2+ntop)+(j+1)*(2+nright)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(3)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(3)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(3)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(3)-1);

        idx = (k+1)*(2+nright)*(2+ntop)+ j   *(2+nright)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(4)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(4)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(4)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(4)-1);

        idx = (k+1)*(2+nright)*(2+ntop)+ j   *(2+nright)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(5)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(5)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(5)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(5)-1);

        idx = (k+1)*(2+nright)*(2+ntop)+(j+1)*(2+nright)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(6)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(6)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(6)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(6)-1);

        idx = (k+1)*(2+nright)*(2+ntop)+(j+1)*(2+nright)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(7)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(7)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(7)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(7)-1);
      }
    mp.front().addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVR, KVT, KVZ, give(C))));
    for (std::size_t v=0; v<variable.size(); v++)
      mp[v+1].addPatch( mp.front()[3].basis().makeGeometry( give(V[v]) ) );

    // Patch #4 (right,bottom):
    volume = domain->searchVolume(4);

    if ((2+nright) * (2+nbottom) * (2+nz) != volume->getNodeNumber())
      GISMO_ERROR("Number of grid points does not match configuration.");
    C.resize((2+nright) * (2+nbottom) * (2+nz), 3);
    for (std::size_t v=0; v<variable.size(); v++)
      V[v].resize((2+nright) * (2+nbottom) * (2+nz), 1);

    elem=0;
    for (auto elemit = volume->element_cbegin(); elemit != volume->element_cend(); elemit++, elem++)
      {
        i = elem%(1+nright);
        j = elem/(1+nright);
        k = elem/((1+nright)*(1+nbottom));

        idx = k   *(2+nright)*(2+nbottom)+ j   *(2+nright)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(0)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(0)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(0)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(0)-1);

        idx = k   *(2+nright)*(2+nbottom)+ j   *(2+nright)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(1)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(1)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(1)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(1)-1);

        idx = k   *(2+nright)*(2+nbottom)+(j+1)*(2+nright)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(2)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(2)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(2)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(2)-1);

        idx = k   *(2+nright)*(2+nbottom)+(j+1)*(2+nright)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(3)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(3)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(3)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(3)-1);

        idx = (k+1)*(2+nright)*(2+nbottom)+ j   *(2+nright)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(4)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(4)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(4)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(4)-1);

        idx = (k+1)*(2+nright)*(2+nbottom)+ j   *(2+nright)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(5)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(5)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(5)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(5)-1);

        idx = (k+1)*(2+nright)*(2+nbottom)+(j+1)*(2+nright)+i;
        C(idx,0) = domain->getNode(elemit->getNodeID(6)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(6)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(6)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(6)-1);

        idx = (k+1)*(2+nright)*(2+nbottom)+(j+1)*(2+nright)+i+1;
        C(idx,0) = domain->getNode(elemit->getNodeID(7)-1).getX();
        C(idx,1) = domain->getNode(elemit->getNodeID(7)-1).getY();
        C(idx,2) = domain->getNode(elemit->getNodeID(7)-1).getZ();
        for (std::size_t v=0; v<variable.size(); v++)
          V[v](idx,0) = variable[v]->getValue(elemit->getNodeID(7)-1);
      }
    mp.front().addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVR, KVB, KVZ, give(C))));
    for (std::size_t v=0; v<variable.size(); v++)
      mp[v+1].addPatch( mp.front()[4].basis().makeGeometry( give(V[v]) ) );

    mp.front().computeTopology();
    mp.front().addAutoBoundaries();

    return mp;
  }

private:
  /// \brief Returns command line with standard arguments set
  gsCmdLine getCmdLine(const std::string& message,
                       const char delimiter = ' ',
                       bool helpAndVersion = true)
  {
    gsCmdLine cmd(message, delimiter, helpAndVersion);

    cmd.addInt("L", "nleft",
               "Number of inner grid points in left part (x-direction).",
               nleft);

    cmd.addInt("M", "nmiddle",
               "Number of inner grid points in middle part (x-direction).",
               nmiddle);

    cmd.addInt("R", "nright",
               "Number of inner grid points in right part (x-direction).",
               nright);

    cmd.addInt("T", "ntop",
               "Number of inner grid points in top part (y-direction).",
               ntop);

    cmd.addInt("B", "nbottom",
               "Number of inner grid points in bottom part (y-direction).",
               nbottom);

    cmd.addInt("Z", "nz",
               "Number of inner grid points in z-direction.",
               nz);

    cmd.addReal("x", "xscale",
                "Scaling factor in x-direction (default x=1.0).",
                xscale);

    cmd.addReal("y", "yscale",
                "Scaling factor in y-direction (default y=1.0).",
                yscale);

    cmd.addReal("z", "zscale",
               "Scaling factor in z-direction (default z=1.0).",
               zscale);

    cmd.addString("", "xleftbottom",
                  "Reparameterization function in left-bottom part (x-direction).",
                  xleftbottom);

    cmd.addString("", "xmiddlebottom",
                  "Reparameterization function in middle-bottom part (x-direction).",
                  xmiddlebottom);

    cmd.addString("", "xrightbottom",
                  "Reparameterization function in right-bottom part (x-direction).",
                  xrightbottom);

    cmd.addString("", "xleftmiddle",
                  "Reparameterization function in left-middle part (x-direction).",
                  xleftmiddle);

    cmd.addString("", "xrightmiddle",
                  "Reparameterization function in right-middle part (x-direction).",
                  xrightmiddle);

    cmd.addString("", "xlefttop",
                  "Reparameterization function in left-top part (x-direction).",
                  xlefttop);

    cmd.addString("", "xmiddletop",
                  "Reparameterization function in middle-top part (x-direction).",
                  xmiddletop);

    cmd.addString("", "xrighttop",
                  "Reparameterization function in right-top part (x-direction).",
                  xrighttop);

    cmd.addString("", "yleftbottomleft",
                  "Reparameterization function in left-bottom part, left boundary (y-direction).",
                  yleftbottomleft);

    cmd.addString("", "yleftbottomright",
                  "Reparameterization function in left-bottom part, right boundary (y-direction).",
                  yleftbottomright);

    cmd.addString("", "ylefttopleft",
                  "Reparameterization function in left-top part, left boundary (y-direction).",
                  ylefttopleft);

    cmd.addString("", "ylefttopright",
                  "Reparameterization function in left-top part, right boundary (y-direction).",
                  ylefttopright);

    cmd.addString("", "yrighttopleft",
                  "Reparameterization function in right-top part, left boundary (y-direction).",
                  yrighttopleft);

    cmd.addString("", "yrighttopright",
                  "Reparameterization function in right-top part, right boundary (y-direction).",
                  yrighttopright);

    cmd.addString("", "yrightbottomright",
                  "Reparameterization function in right-bottom part, right boundary (y-direction).",
                  yrightbottomright);

    cmd.addString("", "yrightbottomleft",
                  "Reparameterization function in right-bottom part, left boundary (y-direction).",
                  yrightbottomleft);

     cmd.addInt( "e", "degreeElevation",
                 "Number of degree elevation steps (default e=0)", numElevate);

     cmd.addInt( "r", "uniformRefine",
                 "Number of uniform h-refinement steps (default r=0)",  numRefine);

     return cmd;
  }

  /// \brief Returns option list with all arguments imported from
  /// file; sets the standard arguments from the imported values (if
  /// available) and the default values otherwise, unless they have
  /// been specified explicitly as command line arguments
  ///
  /// \note The priority of arguments is as follows:
  ///
  /// 1. Command line argument overrides default value and value from parameter file
  /// 2. Value from parameter file overrides default value
  /// 3. Default value overrides invalidated initial value
  gsOptionList getOptionList(const std::string& filename)
  {
    gsOptionList options;

    if (!filename.empty()) {
      gsFileData<> fd(filename);
      fd.getFirst(options);
    }

    if (nleft      == -1) nleft      = options.askInt("nleft",      0);
    if (nmiddle    == -1) nmiddle    = options.askInt("nmiddle",    0);
    if (nright     == -1) nright     = options.askInt("nright",     0);
    if (ntop       == -1) ntop       = options.askInt("ntop",       0);
    if (nbottom    == -1) nbottom    = options.askInt("nbottom",    0);
    if (nz         == -1) nz         = options.askInt("nz",         0);
    if (numElevate == -1) numElevate = options.askInt("numElevate", 0);
    if (numRefine  == -1) numRefine  = options.askInt("numRefine",  0);

    if (xscale   == -1.0) xscale     = options.askReal("xscale",  1.0);
    if (yscale   == -1.0) yscale     = options.askReal("yscale",  1.0);
    if (zscale   == -1.0) zscale     = options.askReal("zscale",  1.0);

    if (xleftbottom      == "INVALID")
      xleftbottom         = options.askString("xleftbottom",   "x");
    if (xleftmiddle      == "INVALID")
      xleftmiddle         = options.askString("xleftmiddle",   "x");
    if (xlefttop         == "INVALID")
      xlefttop            = options.askString("xlefttop",      "x");
    if (xrightbottom     == "INVALID")
      xrightbottom        = options.askString("xrightbottom",  "x");
    if (xrightmiddle     == "INVALID")
      xrightmiddle        = options.askString("xrightmiddle",  "x");
    if (xrighttop        == "INVALID")
      xrighttop           = options.askString("xrighttop",     "x");
    if (xmiddletop       == "INVALID")
      xmiddletop          = options.askString("xmiddletop",    "x");
    if (xmiddlebottom    == "INVALID")
      xmiddlebottom       = options.askString("xmiddlebottom", "x");

    if(yleftbottomleft   == "INVALID")
      yleftbottomleft     = options.askString("yleftbottomleft",   "y");
    if(yleftbottomright  == "INVALID")
      yleftbottomright    = options.askString("yleftbottomright",  "y");
    if(yrightbottomleft  == "INVALID")
      yrightbottomleft    = options.askString("yrightbottomleft",  "y");
    if(yrightbottomright == "INVALID")
      yrightbottomright   = options.askString("yrightbottomright", "y");
    if(ylefttopleft      == "INVALID")
      ylefttopleft        = options.askString("ylefttopleft",      "y");
    if(ylefttopright     == "INVALID")
      ylefttopright       = options.askString("ylefttopright",     "y");
    if(yrighttopleft     == "INVALID")
      yrighttopleft       = options.askString("yrighttopleft",     "y");
    if(yrighttopright    == "INVALID")
      yrighttopright      = options.askString("yrighttopright",    "y");

    if (hash == 0)  hash = options.askInt("hash",  0);

    return options;
  }

  /// \brief Returns option list with all standard arguments set to
  /// the internal values
  gsOptionList setOptionList() const
  {
    gsOptionList options;

    options.addInt("nleft", "Number of inner grid points in left part (x-direction)", nleft);
    options.addInt("nmiddle", "Number of inner grid points in middle part (x-direction)", nmiddle);
    options.addInt("nright", "Number of inner grid points in right part (x-direction)", nright);
    options.addInt("ntop", "Number of inner grid points in top part (y-direction)", ntop);
    options.addInt("nbottom", "Number of inner grid points in bottom part (y-direction)", nbottom);
    options.addInt("nz", "Number of inner grid points in z-direction", nz);
    options.addInt("numElevate", "Number of degree elevation steps", numElevate);
    options.addInt("numRefine", "Number of uniform h-refinement steps", numRefine);

    options.addReal("xscale", "Scaling factor in x-direction", xscale);
    options.addReal("yscale", "Scaling factor in y-direction", yscale);
    options.addReal("zscale", "Scaling factor in z-direction", zscale);

    options.addString("xleftbottom",
                      "Reparameterization function in left-bottom part (x-direction)",
                      xleftbottom);
    options.addString("xmiddlebottom",
                      "Reparameterization function in middle-bottom part (x-direction)",
                      xmiddlebottom);
    options.addString("xrightbottom",
                      "Reparameterization function in right-bottom part (x-direction)",
                      xrightbottom);
    options.addString("xleftmiddle",
                      "Reparameterization function in left-middle part (x-direction)",
                      xleftmiddle);
    options.addString("xrightmiddle",
                      "Reparameterization function in right-middle part (x-direction)",
                      xrightmiddle);
    options.addString("xlefttop",
                      "Reparameterization function in left-top part (x-direction)",
                      xlefttop);
    options.addString("xmiddletop",
                      "Reparameterization function in middle-top part (x-direction)",
                      xmiddletop);
    options.addString("xrighttop",
                      "Reparameterization function in right-top part (x-direction)",
                      xrighttop);
    options.addString("yleftbottomleft",
                      "Reparameterization function in left-bottom part, left boundary (y-direction)",
                      yleftbottomleft);
    options.addString("yleftbottomright",
                      "Reparameterization function in left-bottom part, right boundary (y-direction)",
                      yleftbottomright);
    options.addString("ylefttopleft",
                      "Reparameterization function in left-top part, left boundary (y-direction)",
                      ylefttopleft);
    options.addString("ylefttopright",
                      "Reparameterization function in left-top part, right boundary (y-direction)",
                      ylefttopright);
    options.addString("yrighttopleft",
                      "Reparameterization function in right-top part, left boundary (y-direction)",
                      yrighttopleft);
    options.addString("yrighttopright",
                      "Reparameterization function in right-top part, right boundary (y-direction)",
                      yrighttopright);
    options.addString("yrightbottomright",
                      "Reparameterization function in right-bottom part, right boundary (y-direction)",
                      yrightbottomright);
    options.addString("yrightbottomleft",
                      "Reparameterization function in right-bottom part, left boundary (y-direction)",
                      yrightbottomleft);

    std::size_t hash = getHash();
    if (hash != std::hash<std::string>{}("0"))
      {
        options.addString("hash", "Hash value of generated CFX5 data", util::to_string(hash));
        
        for (auto domit = domain_cbegin(); domit != domain_cend(); domit++)
          {
            hash = (*domit)->getHash();
            options.addString("hash_domain" + util::to_string((*domit)->getDomain()),
                              "Hash value of CFX5 domain "  + (*domit)->getName(), util::to_string(hash));

            for (auto regit = (*domit)->region_cbegin(); regit != (*domit)->region_cend(); regit++)
              {
                hash = (*regit)->getHash();
                options.addString("hash_domain" + util::to_string((*domit)->getDomain())
                                  +   "_region" + util::to_string((*regit)->getRegion()),
                                  "Hash value of CFX5 region "  + (*regit)->getName(), util::to_string(hash));
              }

            for (auto volit = (*domit)->volume_cbegin(); volit != (*domit)->volume_cend(); volit++)
              {
                hash = (*volit)->getHash();
                options.addString("hash_domain" + util::to_string((*domit)->getDomain())
                                  +   "_volume" + util::to_string((*volit)->getVolume()),
                                  "Hash value of CFX5 volume"   + (*volit)->getName(), util::to_string(hash));
              }

            for (auto bdrit = (*domit)->boundary_cbegin(); bdrit != (*domit)->boundary_cend(); bdrit++)
              {
                hash = (*bdrit)->getHash();
                options.addString("hash_domain" + util::to_string((*domit)->getDomain())
                                  +   "_boundary" + util::to_string((*bdrit)->getBoundary()),
                                  "Hash value of CFX5 boundary"   + (*bdrit)->getName(), util::to_string(hash));
              }

            for (auto varit = (*domit)->variable_cbegin(); varit != (*domit)->variable_cend(); varit++)
              {
                hash = (*varit)->getHash();
                options.addString("hash_domain" + util::to_string((*domit)->getDomain())
                                  +   "_variable" + util::to_string((*varit)->getVariable()),
                                  "Hash value of CFX5 variable"   + (*varit)->getName(), util::to_string(hash));
              }

          }
      }

    return options;
  }

  /// \brief Extracts the file name
  static std::string getFileName(const std::string& fullpath)
  {
    size_t i = fullpath.rfind('.', fullpath.length());
    if (i != std::string::npos)
      return (fullpath.substr(0, i));
    else
      return("NOT_FOUND");
  }

  /// \brief Extracts the file extension
  static std::string getFileExt(const std::string& fullpath)
  {
    size_t i = fullpath.rfind('.', fullpath.length());
    if (i != std::string::npos)
      return (fullpath.substr(i+1, fullpath.length() - i));
    else      
      return("NOT_FOUND");
  }

  /// \brief Generates structured multi-patch grid (stored as G+Smo
  /// multi-patch geometry with degree 1) from 2d G+Smo geometry
  gsMultiPatch<> genGridFromGeo2d(const gsMultiPatch<>& geo) const
  {
    // Coefficients
    gsMatrix<> C[5];

    // Patch #0 (left,bottom):
    C[0] = gsMatrix<>(2, (2+nleft) * (2+nbottom));

#pragma omp parallel for collapse(2)
    for (int j=0; j<(2+nbottom); j++)
      for (int i=0; i<(2+nleft); i++)
        {
          C[0](0, (2+nleft)*j + i) = i/real_t(1+nleft);
          C[0](1, (2+nleft)*j + i) = j/real_t(1+nbottom);
        }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xleftbottom     + ") + (y)*(" + xleftmiddle      + ")",
                     "(1-x)*(" + yleftbottomleft + ") + (x)*(" + yleftbottomright + ")",
                     2).eval_into(C[0], C[0]);

    // Patch #1 (left,top)
    C[1] = gsMatrix<>(2, (2+nleft) * (2+ntop));

#pragma omp parallel for collapse(2)
    for (int j=0; j<(2+ntop); j++)
      for (int i=0; i<(2+nleft); i++)
        {
          C[1](0, (2+nleft)*j + i) = i/real_t(1+nleft);
          C[1](1, (2+nleft)*j + i) = j/real_t(1+ntop);
        }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xleftmiddle  + ") + (y)*(" + xlefttop      + ")",
                     "(1-x)*(" + ylefttopleft + ") + (x)*(" + ylefttopright + ")",
                     2).eval_into(C[1], C[1]);

    // Patch #2 (middle,top)
    C[2] = gsMatrix<>(2, (2+nmiddle) * (2+ntop));

#pragma omp parallel for collapse(2)
    for (int j=0; j<(2+ntop); j++)
      for (int i=0; i<(2+nmiddle); i++)
        {
          C[2](0, (2+nmiddle)*j + i) = i/real_t(1+nmiddle);
          C[2](1, (2+nmiddle)*j + i) = j/real_t(1+ntop);
        }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xmiddlebottom + ")   + (y)*(" + xmiddletop    + ")",
                     "(1-x)*(" + ylefttopright   + ") + (x)*(" + yrighttopleft + ")",
                     2).eval_into(C[2], C[2]);

    // Patch #3 (right,top)
    C[3] = gsMatrix<>(2, (2+nright) * (2+ntop));

#pragma omp parallel for collapse(2)
    for (int j=0; j<(2+ntop); j++)
      for (int i=0; i<(2+nright); i++)
        {
          C[3](0, (2+nright)*j + i) = i/real_t(1+nright);
          C[3](1, (2+nright)*j + i) = j/real_t(1+ntop);
        }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xrightmiddle + ")  + (y)*(" + xrighttop      + ")",
                     "(1-x)*(" + yrighttopleft + ") + (x)*(" + yrighttopright + ")",
                     2).eval_into(C[3], C[3]);

    // Patch #4 (right,bottom)
    C[4] = gsMatrix<>(2, (2+nright) * (2+nbottom));

#pragma omp parallel for collapse(2)
    for (int j=0; j<(2+nbottom); j++)
      for (int i=0; i<(2+nright); i++)
        {
          C[4](0, (2+nright)*j + i) = i/real_t(1+nright);
          C[4](1, (2+nright)*j + i) = j/real_t(1+nbottom);
        }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xrightbottom     + ") + (y)*(" + xrightmiddle      + ")",
                     "(1-x)*(" + yrightbottomleft + ") + (x)*(" + yrightbottomright + ")",
                     2).eval_into(C[4], C[4]);

    // Calculate physical coordinates of grid points
    for (int i=0; i<5; i++)
      geo[i].eval_into(C[i], C[i]);

    // Generate uniform knot vectors
    gsKnotVector<> KVL(0, 1, nleft,   2);
    gsKnotVector<> KVM(0, 1, nmiddle, 2);
    gsKnotVector<> KVR(0, 1, nright,  2);
    gsKnotVector<> KVT(0, 1, ntop,    2);
    gsKnotVector<> KVB(0, 1, nbottom, 2);

    // Define type 2D-tensor-B-spline pointer
    typedef typename gsTensorBSpline<2>::uPtr TensorBSpline2Ptr;

    // Generate multi-patch parameterization
    gsMultiPatch<> grid;

    grid.addPatch(TensorBSpline2Ptr(new gsTensorBSpline<2>(KVL, KVB,
                                                           give(C[0].transpose()))));
    grid.addPatch(TensorBSpline2Ptr(new gsTensorBSpline<2>(KVL, KVT,
                                                           give(C[1].transpose()))));
    grid.addPatch(TensorBSpline2Ptr(new gsTensorBSpline<2>(KVM, KVT,
                                                           give(C[2].transpose()))));
    grid.addPatch(TensorBSpline2Ptr(new gsTensorBSpline<2>(KVR, KVT,
                                                           give(C[3].transpose()))));
    grid.addPatch(TensorBSpline2Ptr(new gsTensorBSpline<2>(KVR, KVB,
                                                           give(C[4].transpose()))));
    grid.computeTopology();
    grid.addAutoBoundaries();

    return grid;
  }

  /// \brief Generates structured multi-patch grid (stored as G+Smo
  /// multi-patch geometry with degree 1) from 2d G+Smo geometry
  gsMultiPatch<> genGridFromGeo3d(const gsMultiPatch<>& geo) const
  {
    // Coefficients
    gsMatrix<> C[5];

    // Patch #0 (left,bottom):
    C[0] = gsMatrix<>(3, (2+nz) * (2+nleft) * (1+nbottom));

#pragma omp parallel for collapse(3)
    for (int k=0; k<(2+nz); k++)
      for (int j=0; j<(2+nbottom); j++)
        for (int i=0; i<(2+nleft); i++)
          {
            C[0](0, (2+nleft)*(2+nbottom)*k + (2+nleft)*j + i) = i/real_t(1+nleft);
            C[0](1, (2+nleft)*(2+nbottom)*k + (2+nleft)*j + i) = j/real_t(1+nbottom);
            C[0](2, (2+nleft)*(2+nbottom)*k + (2+nleft)*j + i) = k/real_t(1+nz);
          }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xleftbottom     + ") + (y)*(" + xleftmiddle      + ")",
                     "(1-x)*(" + yleftbottomleft + ") + (x)*(" + yleftbottomright + ")",
                     "z",
                     3).eval_into(C[0], C[0]);

    // Patch #1 (left,top)
    C[1] = gsMatrix<>(3, (2+nz) * (2+nleft) * (2+ntop));

#pragma omp parallel for collapse(3)
    for (int k=0; k<(2+nz); k++)
      for (int j=0; j<(2+ntop); j++)
        for (int i=0; i<(2+nleft); i++)
          {
            C[1](0, (2+nleft)*(2+ntop)*k + (2+nleft)*j + i) = i/real_t(1+nleft);
            C[1](1, (2+nleft)*(2+ntop)*k + (2+nleft)*j + i) = j/real_t(1+ntop);
            C[1](2, (2+nleft)*(2+ntop)*k + (2+nleft)*j + i) = k/real_t(1+nz);
          }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xleftmiddle  + ") + (y)*(" + xlefttop      + ")",
                     "(1-x)*(" + ylefttopleft + ") + (x)*(" + ylefttopright + ")",
                     "z",
                     3).eval_into(C[1], C[1]);

    // Patch #2 (middle,top)
    C[2] = gsMatrix<>(3, (2+nz) * (2+nmiddle) * (2+ntop));

#pragma omp parallel for collapse(3)
    for (int k=0; k<(2+nz); k++)
      for (int j=0; j<(2+ntop); j++)
        for (int i=0; i<(2+nmiddle); i++)
          {
            C[2](0, (2+nmiddle)*(2+ntop)*k + (2+nmiddle)*j + i) = i/real_t(1+nmiddle);
            C[2](1, (2+nmiddle)*(2+ntop)*k + (2+nmiddle)*j + i) = j/real_t(1+ntop);
            C[2](2, (2+nmiddle)*(2+ntop)*k + (2+nmiddle)*j + i) = k/real_t(1+nz);
          }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xmiddlebottom + ")   + (y)*(" + xmiddletop    + ")",
                     "(1-x)*(" + ylefttopright   + ") + (x)*(" + yrighttopleft + ")",
                     "z",
                     3).eval_into(C[2], C[2]);

    // Patch #3 (right,top)
    C[3] = gsMatrix<>(3, (2+nz) * (2+nright) * (2+ntop));

#pragma omp parallel for collapse(2)
    for (int k=0; k<(2+nz); k++)
      for (int j=0; j<(2+ntop); j++)
        for (int i=0; i<(2+nright); i++)
          {
            C[3](0, (2+nright)*(2+ntop)*k + (2+nright)*j + i) = i/real_t(1+nright);
            C[3](1, (2+nright)*(2+ntop)*k + (2+nright)*j + i) = j/real_t(1+ntop);
            C[3](2, (2+nright)*(2+ntop)*k + (2+nright)*j + i) = k/real_t(1+nz);
          }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xrightmiddle + ")  + (y)*(" + xrighttop      + ")",
                     "(1-x)*(" + yrighttopleft + ") + (x)*(" + yrighttopright + ")",
                     "z",
                     3).eval_into(C[3], C[3]);

    // Patch #4 (right,bottom)
    C[4] = gsMatrix<>(3, (2+nz) * (2+nright) * (2+nbottom));

#pragma omp parallel for collapse(2)
    for (int k=0; k<(2+nz); k++)
      for (int j=0; j<(2+nbottom); j++)
        for (int i=0; i<(2+nright); i++)
          {
            C[4](0, (2+nright)*(2+nbottom)*k + (2+nright)*j + i) = i/real_t(1+nright);
            C[4](1, (2+nright)*(2+nbottom)*k + (2+nright)*j + i) = j/real_t(1+nbottom);
            C[4](2, (2+nright)*(2+nbottom)*k + (2+nright)*j + i) = k/real_t(1+nz);
          }

    // Apply reparameterization
    gsFunctionExpr<>("(1-y)*(" + xrightbottom     + ") + (y)*(" + xrightmiddle      + ")",
                     "(1-x)*(" + yrightbottomleft + ") + (x)*(" + yrightbottomright + ")",
                     "z",
                     3).eval_into(C[4], C[4]);

    // Calculate physical coordinates of grid points
    for (int i=0; i<5; i++)
      geo[i].eval_into(C[i], C[i]);

    // Generate uniform knot vectors
    gsKnotVector<> KVZ(0, 1, nz,      2);
    gsKnotVector<> KVL(0, 1, nleft,   2);
    gsKnotVector<> KVM(0, 1, nmiddle, 2);
    gsKnotVector<> KVR(0, 1, nright,  2);
    gsKnotVector<> KVT(0, 1, ntop,    2);
    gsKnotVector<> KVB(0, 1, nbottom, 2);

    // Define type 3D-tensor-B-spline pointer
    typedef typename gsTensorBSpline<3>::uPtr TensorBSpline3Ptr;

    // Generate multi-patch parameterization
    gsMultiPatch<> grid;

    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVL, KVB, KVZ,
                                                           give(C[0].transpose()))));
    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVL, KVT, KVZ,
                                                           give(C[1].transpose()))));
    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVM, KVT, KVZ,
                                                           give(C[2].transpose()))));
    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVR, KVT, KVZ,
                                                           give(C[3].transpose()))));
    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVR, KVB, KVZ,
                                                           give(C[4].transpose()))));
    grid.computeTopology();
    grid.addAutoBoundaries();

    return grid;
  }

  /// \brief Generates CFX5 nodes from 2d grid stored as G+Smo geometry
  void genCFX5NodesFromMultiPatch2d(const gsMultiPatch<> mp, gsCFX5Domain* domain)
  {
    // Compute total number of nodes
    std::size_t nodes2d = ( (nbottom + ntop + 3) * (nleft + nmiddle + nright + 4) -
                            (nbottom + 1)        *  nmiddle
                            );
    std::size_t nodes3d = (2 + nz) * nodes2d;

    // Initialize nodes
    domain->addNodes(nodes3d);

    // Initialize iterator over nodes
    auto nodeit = domain->node_begin();

    // Initialize physical coordinates of grid points
    for (int k=0; k<(2+nz); k++)
      {
        // Patch #0 (left,bottom)
        for (int j=0; j<(1+nbottom); j++)
          for (int i=0; i<(2+nleft); i++)
            {
              nodeit->getX() = xscale * mp[0].coefs()((2+nleft)*j + i, 0);
              nodeit->getY() = yscale * mp[0].coefs()((2+nleft)*j + i, 1);
              nodeit->getZ() = zscale * k/real_t(1+nz);
              nodeit++;
            }

        // Patch #1 (left,top)
        for (int j=0; j<(2+ntop); j++)
          for (int i=0; i<(2+nleft); i++)
            {
              nodeit->getX() = xscale * mp[1].coefs()((2+nleft)*j + i, 0);
              nodeit->getY() = yscale * mp[1].coefs()((2+nleft)*j + i, 1);
              nodeit->getZ() = zscale * k/real_t(1+nz);
              nodeit++;
            }

        // Patch #4 (right,bottom)
        for (int j=0; j<(1+nbottom); j++)
          for (int i=0; i<(2+nright); i++)
            {
              nodeit->getX() = xscale * mp[4].coefs()((2+nright)*j + i, 0);
              nodeit->getY() = yscale * mp[4].coefs()((2+nright)*j + i, 1);
              nodeit->getZ() = zscale * k/real_t(1+nz);
              nodeit++;
            }

        // Patch #3 (right,top)
        for (int j=0; j<(2+ntop); j++)
          for (int i=0; i<(2+nright); i++)
            {
              nodeit->getX() = xscale * mp[3].coefs()((2+nright)*j + i, 0);
              nodeit->getY() = yscale * mp[3].coefs()((2+nright)*j + i, 1);
              nodeit->getZ() = zscale * k/real_t(1+nz);
              nodeit++;
            }

        // Patch #2 (middle,top)
        for (int j=0; j<(2+ntop); j++)
          for (int i=1; i<(1+nmiddle); i++)
            {
              nodeit->getX() = xscale * mp[2].coefs()((2+nmiddle)*j + i, 0);
              nodeit->getY() = yscale * mp[2].coefs()((2+nmiddle)*j + i, 1);
              nodeit->getZ() = zscale * k/real_t(1+nz);
              nodeit++;
            }
      }
  }

  /// \brief Generates CFX5 nodes from 3d grid stored as G+Smo geometry
  void genCFX5NodesFromMultiPatch3d(const gsMultiPatch<> mp, gsCFX5Domain* domain)
  {
    // Compute total number of nodes
    std::size_t nodes2d = ( (nbottom + ntop + 3) * (nleft + nmiddle + nright + 4) -
                            (nbottom + 1)        *  nmiddle
                            );
    std::size_t nodes3d = (2 + nz) * nodes2d;

    // Initialize nodes
    domain->addNodes(nodes3d);

    // Initialize iterator over nodes
    auto nodeit = domain->node_begin();

    // Initialize physical coordinates of grid points
    for (int k=0; k<(2+nz); k++)
      {
        // Patch #0 (left,bottom)
        for (int j=0; j<(1+nbottom); j++)
          for (int i=0; i<(2+nleft); i++)
            {
              nodeit->getX() = xscale * mp[0].coefs()((2+nleft)*(2+nbottom)*k + (2+nleft)*j + i, 0);
              nodeit->getY() = yscale * mp[0].coefs()((2+nleft)*(2+nbottom)*k + (2+nleft)*j + i, 1);
              nodeit->getZ() = zscale * mp[0].coefs()((2+nleft)*(2+nbottom)*k + (2+nleft)*j + i, 2);
              nodeit++;
            }

        // Patch #1 (left,top)
        for (int j=0; j<(2+ntop); j++)
          for (int i=0; i<(2+nleft); i++)
            {
              nodeit->getX() = xscale * mp[1].coefs()((2+nleft)*(2+ntop)*k + (2+nleft)*j + i, 0);
              nodeit->getY() = yscale * mp[1].coefs()((2+nleft)*(2+ntop)*k + (2+nleft)*j + i, 1);
              nodeit->getZ() = zscale * mp[1].coefs()((2+nleft)*(2+ntop)*k + (2+nleft)*j + i, 2);
              nodeit++;
            }

        // Patch #4 (right,bottom)
        for (int j=0; j<(1+nbottom); j++)
          for (int i=0; i<(2+nright); i++)
            {
              nodeit->getX() = xscale * mp[4].coefs()((2+nright)*(2+nbottom)*k + (2+nright)*j + i, 0);
              nodeit->getY() = yscale * mp[4].coefs()((2+nright)*(2+nbottom)*k + (2+nright)*j + i, 1);
              nodeit->getZ() = zscale * mp[4].coefs()((2+nright)*(2+nbottom)*k + (2+nright)*j + i, 2);
              nodeit++;
            }

        // Patch #3 (right,top)
        for (int j=0; j<(2+ntop); j++)
          for (int i=0; i<(2+nright); i++)
            {
              nodeit->getX() = xscale * mp[3].coefs()((2+nright)*(2+ntop)*k + (2+nright)*j + i, 0);
              nodeit->getY() = yscale * mp[3].coefs()((2+nright)*(2+ntop)*k + (2+nright)*j + i, 1);
              nodeit->getZ() = zscale * mp[3].coefs()((2+nright)*(2+ntop)*k + (2+nright)*j + i, 2);
              nodeit++;
            }

        // Patch #2 (middle,top)
        for (int j=0; j<(2+ntop); j++)
          for (int i=1; i<(1+nmiddle); i++)
            {
              nodeit->getX() = xscale * mp[2].coefs()((2+nmiddle)*(2+ntop)*k + (2+nmiddle)*j + i, 0);
              nodeit->getY() = yscale * mp[2].coefs()((2+nmiddle)*(2+ntop)*k + (2+nmiddle)*j + i, 1);
              nodeit->getZ() = zscale * mp[2].coefs()((2+nmiddle)*(2+ntop)*k + (2+nmiddle)*j + i, 2);
              nodeit++;
            }
      }
  }

  /// \brief Generates structured multi-patch CFX5 grid from grid points in 3d
  void genCFX5Elements(gsCFX5Domain* domain)
  {
    // Compute total number of nodes
    std::size_t nodes2d = ( (nbottom + ntop + 3) * (nleft + nmiddle + nright + 4) -
                            (nbottom + 1)        *  nmiddle
                            );

    // Compute total number of elements
    std::size_t elements = (1 + nz) * ( (nbottom + ntop + 2) * (nleft + nmiddle + nright + 3) -
                                        (nbottom + 1)        * (nmiddle + 1)
                                        );

    // Initialize elements
    domain->addElements(elements);

    // Initialize iterator over nodes and elements
    auto elemit = domain->element_begin();
    int element = 0;

    // Initialize patch #0 (left,bottom)
    gsCFX5Volume* vol0 = domain->addVolume(domain);
    vol0->setVolume(0);
    vol0->setName("Fluid0");
    vol0->addElementIDs( (1+nz)*(1+nbottom)*(1+nleft) );
    auto vol0it = vol0->elementid_begin();

    gsCFX5Region* reg0north = domain->addRegion(domain);
    gsCFX5Region* reg0south = domain->addRegion(domain);
    gsCFX5Region* reg0east  = domain->addRegion(domain);
    gsCFX5Region* reg0west  = domain->addRegion(domain);
    gsCFX5Region* reg0front = domain->addRegion(domain);
    gsCFX5Region* reg0back  = domain->addRegion(domain);
    reg0north->setRegion(0);
    reg0south->setRegion(1);
    reg0east-> setRegion(2);
    reg0west-> setRegion(3);
    reg0front->setRegion(4);
    reg0back-> setRegion(5);
    reg0north->setName("Fluid0_North");
    reg0south->setName("Fluid0_South");
    reg0east-> setName("Fluid0_East");
    reg0west-> setName("Fluid0_West");
    reg0front->setName("Fluid0_Front");
    reg0back-> setName("Fluid0_Back");
    reg0north->addFaceIDs( (1+nz)*(1+nleft) );
    reg0south->addFaceIDs( (1+nz)*(1+nleft) );
    reg0east-> addFaceIDs( (1+nz)*(1+nbottom) );
    reg0west-> addFaceIDs( (1+nz)*(1+nbottom) );
    reg0front->addFaceIDs( (1+nleft)*(1+nbottom) );
    reg0back-> addFaceIDs( (1+nleft)*(1+nbottom) );
    auto reg0northit = reg0north->faceid_begin();
    auto reg0southit = reg0south->faceid_begin();
    auto reg0eastit  = reg0east-> faceid_begin();
    auto reg0westit  = reg0west-> faceid_begin();
    auto reg0frontit = reg0front->faceid_begin();
    auto reg0backit  = reg0back-> faceid_begin();

    // Initialize patch #1 (left,top)
    gsCFX5Volume* vol1 = domain->addVolume(domain);
    vol1->setVolume(1);
    vol1->setName("Fluid1");
    vol1->addElementIDs( (1+nz)*(1+ntop)*(1+nleft) );
    auto vol1it = vol1->elementid_begin();

    gsCFX5Region* reg1north = domain->addRegion(domain);
    gsCFX5Region* reg1south = domain->addRegion(domain);
    gsCFX5Region* reg1east  = domain->addRegion(domain);
    gsCFX5Region* reg1west  = domain->addRegion(domain);
    gsCFX5Region* reg1front = domain->addRegion(domain);
    gsCFX5Region* reg1back  = domain->addRegion(domain);
    reg1north->setRegion(6);
    reg1south->setRegion(7);
    reg1east-> setRegion(8);
    reg1west-> setRegion(9);
    reg1front->setRegion(10);
    reg1back-> setRegion(11);
    reg1north->setName("Fluid1_North");
    reg1south->setName("Fluid1_South");
    reg1east-> setName("Fluid1_East");
    reg1west-> setName("Fluid1_West");
    reg1front->setName("Fluid1_Front");
    reg1back-> setName("Fluid1_Back");
    reg1north->addFaceIDs( (1+nz)*(1+nleft) );
    reg1south->addFaceIDs( (1+nz)*(1+nleft) );
    reg1east-> addFaceIDs( (1+nz)*(1+ntop) );
    reg1west-> addFaceIDs( (1+nz)*(1+ntop) );
    reg1front->addFaceIDs( (1+nleft)*(1+ntop) );
    reg1back-> addFaceIDs( (1+nleft)*(1+ntop) );
    auto reg1northit = reg1north->faceid_begin();
    auto reg1southit = reg1south->faceid_begin();
    auto reg1eastit  = reg1east-> faceid_begin();
    auto reg1westit  = reg1west-> faceid_begin();
    auto reg1frontit = reg1front->faceid_begin();
    auto reg1backit  = reg1back-> faceid_begin();

    // Initialize patch #2 (middle,top)
    gsCFX5Volume* vol2 = domain->addVolume(domain);
    vol2->setVolume(2);
    vol2->setName("Fluid2");
    vol2->addElementIDs( (1+nz)*(1+ntop)*(1+nmiddle) );
    auto vol2it = vol2->elementid_begin();

    gsCFX5Region* reg2north = domain->addRegion(domain);
    gsCFX5Region* reg2south = domain->addRegion(domain);
    gsCFX5Region* reg2east  = domain->addRegion(domain);
    gsCFX5Region* reg2west  = domain->addRegion(domain);
    gsCFX5Region* reg2front = domain->addRegion(domain);
    gsCFX5Region* reg2back  = domain->addRegion(domain);
    reg2north->setRegion(12);
    reg2south->setRegion(13);
    reg2east-> setRegion(14);
    reg2west-> setRegion(15);
    reg2front->setRegion(16);
    reg2back-> setRegion(17);
    reg2north->setName("Fluid2_North");
    reg2south->setName("Fluid2_South");
    reg2east-> setName("Fluid2_East");
    reg2west-> setName("Fluid2_West");
    reg2front->setName("Fluid2_Front");
    reg2back-> setName("Fluid2_Back");
    reg2north->addFaceIDs( (1+nz)*(1+nmiddle) );
    reg2south->addFaceIDs( (1+nz)*(1+nmiddle) );
    reg2east-> addFaceIDs( (1+nz)*(1+ntop) );
    reg2west-> addFaceIDs( (1+nz)*(1+ntop) );
    reg2front->addFaceIDs( (1+nmiddle)*(1+ntop) );
    reg2back-> addFaceIDs( (1+nmiddle)*(1+ntop) );
    auto reg2northit = reg2north->faceid_begin();
    auto reg2southit = reg2south->faceid_begin();
    auto reg2eastit  = reg2east-> faceid_begin();
    auto reg2westit  = reg2west-> faceid_begin();
    auto reg2frontit = reg2front->faceid_begin();
    auto reg2backit  = reg2back-> faceid_begin();

    // Initialize patch #3 (right,top)
    gsCFX5Volume* vol3 = domain->addVolume(domain);
    vol3->setVolume(3);
    vol3->setName("Fluid3");
    vol3->addElementIDs( (1+nz)*(1+ntop)*(1+nright) );
    auto vol3it = vol3->elementid_begin();

    gsCFX5Region* reg3north = domain->addRegion(domain);
    gsCFX5Region* reg3south = domain->addRegion(domain);
    gsCFX5Region* reg3east  = domain->addRegion(domain);
    gsCFX5Region* reg3west  = domain->addRegion(domain);
    gsCFX5Region* reg3front = domain->addRegion(domain);
    gsCFX5Region* reg3back  = domain->addRegion(domain);
    reg3north->setRegion(18);
    reg3south->setRegion(19);
    reg3east-> setRegion(20);
    reg3west-> setRegion(21);
    reg3front->setRegion(22);
    reg3back-> setRegion(23);
    reg3north->setName("Fluid3_North");
    reg3south->setName("Fluid3_South");
    reg3east-> setName("Fluid3_East");
    reg3west-> setName("Fluid3_West");
    reg3front->setName("Fluid3_Front");
    reg3back-> setName("Fluid3_Back");
    reg3north->addFaceIDs( (1+nz)*(1+nright) );
    reg3south->addFaceIDs( (1+nz)*(1+nright) );
    reg3east-> addFaceIDs( (1+nz)*(1+ntop) );
    reg3west-> addFaceIDs( (1+nz)*(1+ntop) );
    reg3front->addFaceIDs( (1+nright)*(1+ntop) );
    reg3back-> addFaceIDs( (1+nright)*(1+ntop) );
    auto reg3northit = reg3north->faceid_begin();
    auto reg3southit = reg3south->faceid_begin();
    auto reg3eastit  = reg3east-> faceid_begin();
    auto reg3westit  = reg3west-> faceid_begin();
    auto reg3frontit = reg3front->faceid_begin();
    auto reg3backit  = reg3back-> faceid_begin();

    // Initialize patch #4 (right,bottom)
    gsCFX5Volume* vol4 = domain->addVolume(domain);
    vol4->setVolume(4);
    vol4->setName("Fluid4");
    vol4->addElementIDs( (1+nz)*(1+nbottom)*(1+nright) );
    auto vol4it = vol4->elementid_begin();

    gsCFX5Region* reg4north = domain->addRegion(domain);
    gsCFX5Region* reg4south = domain->addRegion(domain);
    gsCFX5Region* reg4east  = domain->addRegion(domain);
    gsCFX5Region* reg4west  = domain->addRegion(domain);
    gsCFX5Region* reg4front = domain->addRegion(domain);
    gsCFX5Region* reg4back  = domain->addRegion(domain);
    reg4north->setRegion(24);
    reg4south->setRegion(25);
    reg4east-> setRegion(26);
    reg4west-> setRegion(27);
    reg4front->setRegion(28);
    reg4back-> setRegion(29);
    reg4north->setName("Fluid4_North");
    reg4south->setName("Fluid4_South");
    reg4east-> setName("Fluid4_East");
    reg4west-> setName("Fluid4_West");
    reg4front->setName("Fluid4_Front");
    reg4back-> setName("Fluid4_Back");
    reg4north->addFaceIDs( (1+nz)*(1+nright) );
    reg4south->addFaceIDs( (1+nz)*(1+nright) );
    reg4east-> addFaceIDs( (1+nz)*(1+nbottom) );
    reg4west-> addFaceIDs( (1+nz)*(1+nbottom) );
    reg4front->addFaceIDs( (1+nright)*(1+nbottom) );
    reg4back-> addFaceIDs( (1+nright)*(1+nbottom) );
    auto reg4northit = reg4north->faceid_begin();
    auto reg4southit = reg4south->faceid_begin();
    auto reg4eastit  = reg4east-> faceid_begin();
    auto reg4westit  = reg4west-> faceid_begin();
    auto reg4frontit = reg4front->faceid_begin();
    auto reg4backit  = reg4back-> faceid_begin();

    // Initialize hexahedral elements
    for (int k=0; k<(1+nz); k++)
      {
        // Patch #0 (left,bottom)
        for (int j=0; j<(1+nbottom); j++)
          for (int i=0; i<(1+nleft); i++)
            {
              elemit->setElementType(gsCFX5ElementType::Hexahedral);
              elemit->getNodeID(0) = 1 + nodes2d* k    + (2+nleft)* j    + i;
              elemit->getNodeID(1) = 1 + nodes2d* k    + (2+nleft)* j    + i+1;
              elemit->getNodeID(2) = 1 + nodes2d* k    + (2+nleft)*(j+1) + i;
              elemit->getNodeID(3) = 1 + nodes2d* k    + (2+nleft)*(j+1) + i+1;

              elemit->getNodeID(4) = 1 + nodes2d*(k+1) + (2+nleft)* j    + i;
              elemit->getNodeID(5) = 1 + nodes2d*(k+1) + (2+nleft)* j    + i+1;
              elemit->getNodeID(6) = 1 + nodes2d*(k+1) + (2+nleft)*(j+1) + i;
              elemit->getNodeID(7) = 1 + nodes2d*(k+1) + (2+nleft)*(j+1) + i+1;

              if (!isValidElement(*elemit, domain))
                GISMO_ERROR("Invalid element: " << (*elemit) );

              if (i==0)       { *reg0westit  = (((1+element) << 3) | 1); reg0westit++;  }
              if (i==nleft)   { *reg0eastit  = (((1+element) << 3) | 2); reg0eastit++;  }
              if (j==0)       { *reg0southit = (((1+element) << 3) | 3); reg0southit++; }
              if (j==nbottom) { *reg0northit = (((1+element) << 3) | 4); reg0northit++; }
              if (k==0)       { *reg0backit  = (((1+element) << 3) | 5); reg0backit++;  }
              if (k==nz)      { *reg0frontit = (((1+element) << 3) | 6); reg0frontit++; }

              *vol0it = 1 + element++; vol0it++; elemit++;
            }

        // Patch #1 (left,top)
        for (int j=0; j<(1+ntop); j++)
          for (int i=0; i<(1+nleft); i++)
            {
              elemit->setElementType(gsCFX5ElementType::Hexahedral);
              elemit->getNodeID(0) = 1 + (2+nleft)*(1+nbottom)
                + nodes2d* k    + (2+nleft)* j    + i;
              elemit->getNodeID(1) = 1 + (2+nleft)*(1+nbottom)
                + nodes2d* k    + (2+nleft)* j    + i+1;
              elemit->getNodeID(2) = 1 + (2+nleft)*(1+nbottom)
                + nodes2d* k    + (2+nleft)*(j+1) + i;
              elemit->getNodeID(3) = 1 + (2+nleft)*(1+nbottom) + nodes2d* k
                + (2+nleft)*(j+1) + i+1;

              elemit->getNodeID(4) = 1 + (2+nleft)*(1+nbottom)
                + nodes2d*(k+1) + (2+nleft)* j    + i;
              elemit->getNodeID(5) = 1 + (2+nleft)*(1+nbottom)
                + nodes2d*(k+1) + (2+nleft)* j    + i+1;
              elemit->getNodeID(6) = 1 + (2+nleft)*(1+nbottom)
                + nodes2d*(k+1) + (2+nleft)*(j+1) + i;
              elemit->getNodeID(7) = 1 + (2+nleft)*(1+nbottom)
                + nodes2d*(k+1) + (2+nleft)*(j+1) + i+1;

              if (!isValidElement(*elemit, domain))
                GISMO_ERROR("Invalid element: " << (*elemit) );

              if (i==0)     { *reg1westit  = (((1+element) << 3) | 1); reg1westit++;  }
              if (i==nleft) { *reg1eastit  = (((1+element) << 3) | 2); reg1eastit++;  }
              if (j==0)     { *reg1southit = (((1+element) << 3) | 3); reg1southit++; }
              if (j==ntop)  { *reg1northit = (((1+element) << 3) | 4); reg1northit++; }
              if (k==0)     { *reg1backit  = (((1+element) << 3) | 5); reg1backit++;  }
              if (k==nz)    { *reg1frontit = (((1+element) << 3) | 6); reg1frontit++; }

              *vol1it = 1 + element++; vol1it++; elemit++;
            }

        // Patch #4 (right,bottom)
        for (int j=0; j<(1+nbottom); j++)
          for (int i=0; i<(1+nright); i++)
            {
              elemit->setElementType(gsCFX5ElementType::Hexahedral);
              elemit->getNodeID(0) = 1 + (2+nleft)*(3+nbottom+ntop)
                + nodes2d* k    + (2+nright)* j    + i;
              elemit->getNodeID(1) = 1 + (2+nleft)*(3+nbottom+ntop)
                + nodes2d* k    + (2+nright)* j    + i+1;
              elemit->getNodeID(2) = 1 + (2+nleft)*(3+nbottom+ntop)
                + nodes2d* k    + (2+nright)*(j+1) + i;
              elemit->getNodeID(3) = 1 + (2+nleft)*(3+nbottom+ntop)
                + nodes2d* k    + (2+nright)*(j+1) + i+1;

              elemit->getNodeID(4) = 1 + (2+nleft)*(3+nbottom+ntop)
                + nodes2d*(k+1) + (2+nright)* j    + i;
              elemit->getNodeID(5) = 1 + (2+nleft)*(3+nbottom+ntop)
                + nodes2d*(k+1) + (2+nright)* j    + i+1;
              elemit->getNodeID(6) = 1 + (2+nleft)*(3+nbottom+ntop)
                + nodes2d*(k+1) + (2+nright)*(j+1) + i;
              elemit->getNodeID(7) = 1 + (2+nleft)*(3+nbottom+ntop)
                + nodes2d*(k+1) + (2+nright)*(j+1) + i+1;

              if (!isValidElement(*elemit, domain))
                GISMO_ERROR("Invalid element: " << (*elemit) );
              
              if (i==0)       { *reg4westit  = (((1+element) << 3) | 1); reg4westit++;  }
              if (i==nright)  { *reg4eastit  = (((1+element) << 3) | 2); reg4eastit++;  }
              if (j==0)       { *reg4southit = (((1+element) << 3) | 3); reg4southit++; }
              if (j==nbottom) { *reg4northit = (((1+element) << 3) | 4); reg4northit++; }
              if (k==0)       { *reg4backit  = (((1+element) << 3) | 5); reg4backit++;  }
              if (k==nz)      { *reg4frontit = (((1+element) << 3) | 6); reg4frontit++; }

              *vol4it = 1 + element++; vol4it++; elemit++;
            }

        // Patch #3 (right,top)
        for (int j=0; j<(1+ntop); j++)
          for (int i=0; i<(1+nright); i++)
            {
              elemit->setElementType(gsCFX5ElementType::Hexahedral);
              elemit->getNodeID(0) = 1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                + nodes2d* k    + (2+nright)* j    + i;
              elemit->getNodeID(1) = 1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                + nodes2d* k    + (2+nright)* j    + i+1;
              elemit->getNodeID(2) = 1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                + nodes2d* k    + (2+nright)*(j+1) + i;
              elemit->getNodeID(3) = 1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                + nodes2d* k    + (2+nright)*(j+1) + i+1;

              elemit->getNodeID(4) = 1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                + nodes2d*(k+1) + (2+nright)* j    + i;
              elemit->getNodeID(5) = 1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                + nodes2d*(k+1) + (2+nright)* j    + i+1;
              elemit->getNodeID(6) = 1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                + nodes2d*(k+1) + (2+nright)*(j+1) + i;
              elemit->getNodeID(7) = 1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                + nodes2d*(k+1) + (2+nright)*(j+1) + i+1;
              
              if (!isValidElement(*elemit, domain))
                GISMO_ERROR("Invalid element: " << (*elemit) );

              if (i==0)      { *reg3westit  = (((1+element) << 3) | 1); reg3westit++;  }
              if (i==nright) { *reg3eastit  = (((1+element) << 3) | 2); reg3eastit++;  }
              if (j==0)      { *reg3southit = (((1+element) << 3) | 3); reg3southit++; }
              if (j==ntop)   { *reg3northit = (((1+element) << 3) | 4); reg3northit++; }
              if (k==0)      { *reg3backit  = (((1+element) << 3) | 5); reg3backit++;  }
              if (k==nz)     { *reg3frontit = (((1+element) << 3) | 6); reg3frontit++; }

              *vol3it = 1 + element++; vol3it++; elemit++;
            }

        // Patch #2 (middle,top)
        for (int j=0; j<(1+ntop); j++)
          for (int i=0; i<(1+nmiddle); i++)
            {
              elemit->setElementType(gsCFX5ElementType::Hexahedral);

              elemit->getNodeID(0) = (i==0
                                      ?
                                      1 + (2+nleft)*(1+nbottom)
                                      + nodes2d* k    + (2+nleft)* j   + nleft+1
                                      :
                                      1 + (4+nleft+nright)*(3+nbottom+ntop)
                                      + nodes2d* k    + (nmiddle)* j + i-1);
              elemit->getNodeID(1) = (i==nmiddle
                                      ?
                                      1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                                      + nodes2d* k    + (2+nright)* j
                                      :
                                      1 + (4+nleft+nright)*(3+nbottom+ntop)
                                      + nodes2d* k    + (nmiddle)* j + i);
              elemit->getNodeID(2) = (i==0
                                      ?
                                      1 + (2+nleft)*(1+nbottom) + nodes2d* k
                                      + (2+nleft)*(j+1) + nleft+1
                                      :
                                      1 + (4+nleft+nright)*(3+nbottom+ntop)
                                      + nodes2d* k    + (nmiddle)*(j+1) + i-1);
              elemit->getNodeID(3) = (i==nmiddle
                                      ?
                                      1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                                      + nodes2d* k    + (2+nright)*(j+1)
                                      :
                                      1 + (4+nleft+nright)*(3+nbottom+ntop)
                                      + nodes2d* k    + (nmiddle)*(j+1) + i);

              elemit->getNodeID(4) = (i==0
                                      ?
                                      1 + (2+nleft)*(1+nbottom)
                                      + nodes2d*(k+1) + (2+nleft)* j   + nleft+1
                                      :
                                      1 + (4+nleft+nright)*(3+nbottom+ntop)
                                      + nodes2d*(k+1) + (nmiddle)* j + i-1);
              elemit->getNodeID(5) = (i==nmiddle
                                      ?
                                      1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                                      + nodes2d*(k+1) + (2+nright)* j
                                      :
                                      1 + (4+nleft+nright)*(3+nbottom+ntop)
                                      + nodes2d*(k+1) + (nmiddle)* j + i);
              elemit->getNodeID(6) = (i==0
                                      ?
                                      1 + (2+nleft)*(1+nbottom)
                                      + nodes2d*(k+1) + (2+nleft)*(j+1) + nleft+1
                                      :
                                      1 + (4+nleft+nright)*(3+nbottom+ntop)
                                      + nodes2d*(k+1) + (nmiddle)*(j+1) + i-1);
              elemit->getNodeID(7) = (i==nmiddle
                                      ?
                                      1 + (2+nleft)*(3+nbottom+ntop) + (2+nright)*(1+nbottom)
                                      + nodes2d*(k+1) + (2+nright)*(j+1)
                                      :
                                      1 + (4+nleft+nright)*(3+nbottom+ntop)
                                      + nodes2d*(k+1) + (nmiddle)*(j+1) + i);

              if (!isValidElement(*elemit, domain))
                GISMO_ERROR("Invalid element: " << (*elemit) );

              if (i==0)       { *reg2westit  = (((1+element) << 3) | 1); reg2westit++;  }
              if (i==nmiddle) { *reg2eastit  = (((1+element) << 3) | 2); reg2eastit++;  }
              if (j==0)       { *reg2southit = (((1+element) << 3) | 3); reg2southit++; }
              if (j==ntop)    { *reg2northit = (((1+element) << 3) | 4); reg2northit++; }
              if (k==0)       { *reg2backit  = (((1+element) << 3) | 5); reg2backit++;  }
              if (k==nz)      { *reg2frontit = (((1+element) << 3) | 6); reg2frontit++; }

              *vol2it = 1 + element++; vol2it++; elemit++;
            }
      }

    // Update internal counters
    domain->calcElementTypes();
  }

  /// \brief Generates G+Smo multi-patch of degree 1 from internally
  /// stored structured multi-patch CFX5 mesh and results
  gsMultiPatch<> genMultiPatchFromCFX5()
  {
    // Coefficients
    gsMatrix<> C[5];

    // Generate uniform knot vectors
    gsKnotVector<> KVZ(0, 1, nz,      2);
    gsKnotVector<> KVL(0, 1, nleft,   2);
    gsKnotVector<> KVM(0, 1, nmiddle, 2);
    gsKnotVector<> KVR(0, 1, nright,  2);
    gsKnotVector<> KVT(0, 1, ntop,    2);
    gsKnotVector<> KVB(0, 1, nbottom, 2);

    // Define type 3D-tensor-B-spline pointer
    typedef typename gsTensorBSpline<3>::uPtr TensorBSpline3Ptr;

    // Generate multi-patch parameterization
    gsMultiPatch<> grid;

    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVL, KVB, KVZ,
                                                           give(C[0].transpose()))));
    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVL, KVT, KVZ,
                                                           give(C[1].transpose()))));
    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVM, KVT, KVZ,
                                                           give(C[2].transpose()))));
    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVR, KVT, KVZ,
                                                           give(C[3].transpose()))));
    grid.addPatch(TensorBSpline3Ptr(new gsTensorBSpline<3>(KVR, KVB, KVZ,
                                                           give(C[4].transpose()))));
    grid.computeTopology();
    grid.addAutoBoundaries();

    return grid;
  }

private:
  /// \brief Time measurement
  gsStopwatch timer;

  /// \brief Number of inner points in the different zones
  int nleft=-1, nmiddle=-1, nright=-1, ntop=-1, nbottom=-1, nz=-1;

  /// \brief Scaling factors
  real_t xscale=-1.0, yscale=-1.0, zscale=-1.0;

  /// \brief Reparameterization functions
  std::string
    xleftbottom       = "INVALID",
    xleftmiddle       = "INVALID",
    xlefttop          = "INVALID",
    xrightbottom      = "INVALID",
    xrightmiddle      = "INVALID",
    xrighttop         = "INVALID",
    xmiddletop        = "INVALID",
    xmiddlebottom     = "INVALID",
    yleftbottomleft   = "INVALID",
    yleftbottomright  = "INVALID",
    yrightbottomleft  = "INVALID",
    yrightbottomright = "INVALID",
    ylefttopleft      = "INVALID",
    ylefttopright     = "INVALID",
    yrighttopleft     = "INVALID",
    yrighttopright    = "INVALID";

  /// \brief Number of degree elevation steps
  int numElevate = -1;

  /// \brief Number of uniform refinement steps
  int numRefine = -1;

  /// Hash value of generated CFX5 data
  std::size_t hash = 0;
};

} // namespace gismo

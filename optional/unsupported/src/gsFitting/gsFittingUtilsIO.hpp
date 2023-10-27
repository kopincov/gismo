/** @file gsFittingIOUtils.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingUtilsIO.h>
#include <gsFitting/gsFittingRoutines.hpp>

namespace gismo
{

template<class T>
int getPoints(std::vector< gsMatrix<T> >& params,
              std::vector< gsMatrix<T> >& points,
              const std::string& filename)
{
    gsFileData<> fd(filename);
    int size = fd.count< gsMatrix<> >()/2;
    for(int i = 0;i < size;i++)
    {
        const int id_par = 2 * i;
        const int id_pts = 2 * i + 1;
        gsMatrix<>* par = fd.getId< gsMatrix<> >(id_par).release();
        gsMatrix<T>* pts = fd.getId< gsMatrix<T> >(id_pts).release();
        params.push_back(*par);
        points.push_back(*pts);
        if(par == NULL || pts == NULL)
            return i;
    }
    return size;
}


template<class T>
int readInput(gsFittingParam<T> &fitting_param)
{
    gsTemplateTargetFitting<T> data(fitting_param);
    gsMultiPatch<T>* topology = NULL;
    gsMultiBasis<T>* basis = NULL;
    bool unique_patch;

    /// There is some file indicating the topology. We are in the case of a multipatch fitting
    if(! fitting_param.topo_file.empty())
    {
        gsFileData<> file_topo
            (fitting_param.topo_file);
        if(file_topo.has<gsMultiPatch<T> >())
        {
            topology = file_topo.getAnyFirst
                < gsMultiPatch<T> >().release();
        }
        else
        {
            GISMO_ENSURE(file_topo.has<gsMultiBasis<T> >(),
                         "Error: the topology file does not contain any multipatch or any multibasis");
            basis = file_topo.getAnyFirst
                < gsMultiBasis<T> >().release();

        }
    }
    unique_patch = ( topology == NULL && basis == NULL);

    gsFileData<> filedata(fitting_param.input_file);
    if ( filedata.has< gsMultiPatch<T> >() )
    {

        gsFileData<> file_template
            (fitting_param.template_file);
        gsMultiPatch<T>* geom = filedata.
            getAnyFirst<gsMultiPatch<T> >().release();
        gsMultiPatch<T>* temp = file_template.
            getAnyFirst<gsMultiPatch<T> >().release();
        GISMO_ASSERT(geom != NULL && temp != NULL,
                     "Error: cannot read one of the input files");
        /*   fittingTargetDyn<gsMultiPatch<T>, T>
             (fitting_param, *geom, *temp, topology, basis);  */
        data = gsTemplateTargetFitting<T>(*geom, *temp, unique_patch,
                                          fitting_param);

    }

    else if ( filedata.has< gsGeometry<T> >() )
    {
        gsFileData<> file_template
            (fitting_param.template_file);
        gsGeometry<T>* geom = filedata.
            getAnyFirst<gsGeometry<T> >().release();
        gsGeometry<T>* temp = file_template.
            getAnyFirst<gsGeometry<T> >().release();
        GISMO_ASSERT(geom != NULL && temp != NULL,
                     "Error: cannot read one of the input files");
        /*  fittingTargetDyn<gsGeometry<T>, T>
            (fitting_param, *geom, *temp, topology, basis);  */

        data = gsTemplateTargetFitting<T>(*geom, *temp, fitting_param);
    }

    if ( filedata.has< gsSolid<T> >() )
    {
        gsFileData<> file_template
            (fitting_param.template_file);
        gsSolid<T>* geom = filedata.getFirst< gsSolid<T> >().release();
        gsSolid<T>* temp = file_template.
            getAnyFirst<gsSolid<T> >().release();
        GISMO_ASSERT(geom != NULL && temp != NULL,
                     "Error: cannot read one of the input files");
        /*  fittingSolid<T>(temp, geom, fitting_param);  */

        data = gsTemplateTargetFitting<T>(*geom, *temp, fitting_param);
    }
/* TODO?
    if ( filedata.has< gsTrimSurface<T> >() )
    {
        gsTrimSurface<T>* bb
            = filedata.getFirst< gsTrimSurface<T> >().release();

            }*/
    else if ( filedata.has< gsMatrix<T> >() )
    {
        std::vector< gsMatrix<T> > param;
        std::vector< gsMatrix<T> > pts;
        getPoints<T>(param, pts, fitting_param.input_file);
        GISMO_ASSERT(param.size() == pts.size() && param.size() > 0,
                     "Error while reading the points");
        /*   fittingPointsDyn<T>(fitting_param, pts, param,
             topology, basis);  */

        data = gsTemplateTargetFitting<T>(pts, param, unique_patch,
                                          fitting_param);
    }
    fittingTrimming<T>(topology, basis, data);
    return 0;
}



template<class GenGeom, class T>
void exportGeometry(gsFunctionSet<T>* map, std::string output,
                    bool exp_xml, int NparameterLines,
                    int dir_slices, int nSlices)
{
    int nPoints = 60;
    GISMO_ENSURE(map != NULL, "The mapping to be exported is NULL");
    GenGeom* _map = static_cast<GenGeom*>(map);
    gsInfo << "Created: " << *_map << "\n\n";

    gsInfo << "Saving: " << output << "\n";
    gsWrite(*_map, output);
    std::string parameter_lines_filename
        = output + "_parameter_lines";
    writeParameterLines(*_map, parameter_lines_filename,
                        nPoints, NparameterLines);
    if(nSlices > 0)
        exportSlices(*_map, output, dir_slices, nSlices, nPoints);
    if(exp_xml)
    {
        std::string name = output;
        gsFileData<> fd(name);
        fd << *_map;
    }
}

template<class T>
void exportSlices(gsGeometry<T>& map, std::string output,
                  int dir_slices, int nSlices, int nPoints)
{
    int dim = map.parDim();
    if(dim == 3)
    {
        gsMatrix<T> bb = map.basis().support();
        T init = bb(dir_slices, 0);
        T del = ( bb(dir_slices, 1) - init ) / (nSlices - 1);
        for(int i = 0;i < nSlices;i++)
        {
            std::string slice_filename
                = output + "_slice_" + util::to_string(i);
            gsGeometrySlice<T> slice =
                map.getIsoParametricSlice(dir_slices, i * del);
            gsWriteParaview(slice, slice_filename, nPoints);
        }
    }
}

template<class T>
void exportSlices(gsMultiPatch<T>& map, std::string output,
                  int dir_slices, int nSlices, int nPoints)
{
    unsigned nPatches =map.nPatches();
    for(unsigned patch = 0;patch < nPatches;patch++)
    {
        std::string patch_path
            = output + "_patch_" + util::to_string(patch);
        exportSlices(map.patch(patch), patch_path,
                     dir_slices, nSlices, nPoints);
    }
}


/// taken from motor/jku/gsMotorIOUtils
void writeKnots(const gsGeometry<>& map,
                const std::string& outputPrefix)
{
    gsMesh<> mesh;
    makeMesh<>(map.basis(), mesh, 3);
    std::string out = outputPrefix + "_Knot_config";
    gsWriteParaview(mesh, out);

    map.evaluateMesh(mesh);
    out = outputPrefix + "_Knot_config_physical";
    gsWriteParaview(mesh, out);
}

void writeKnots(const gsGeometry<>& map,
                const int id,
                const std::string& outputPrefix)
{
    gsMesh<> mesh;
    makeMesh<>(map.basis(), mesh, 3);
    std::string out = outputPrefix + "_Knot_config_" + util::to_string(id);
    gsWriteParaview(mesh, out);

    map.evaluateMesh(mesh);
    out = outputPrefix + "_Knot_config_physical_" + util::to_string(id);
    gsWriteParaview(mesh, out);
}


} /// namespace gismo

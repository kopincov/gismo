/** @file gsFittingUtilsSolids.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingUtilsSolids.h>
#include <gsFitting/gsFittingConstr.hpp>

namespace gismo
{

template<class T> gsFunctionSet<T>*
fitSolidFixedParam(gsSolid<T>* temp, gsSolid<T>* geom,
                   gsFittingParam<T>& fitting_param,
                   int nPointsSample)
{
    typedef typename gsSolidElement<T>
        ::gsSolidHalfFaceHandle gsSolidHalfFaceHandle;

    std::vector<gsSolidHalfFaceHandle> faces_temp =
        temp->face;
    std::vector<gsSolidHalfFaceHandle> faces_geom =
        geom->face;
    GISMO_ENSURE(faces_geom.size() == faces_temp.size(),
                 "ERROR: Diffent number of faces");

    unsigned s = faces_geom.size();
    gsFunctionSet<T>* _result = NULL;
    gsGeometry<T>* result = NULL;

    std::vector<std::vector<gsMatrix<T> > > points_loc;
    std::vector<std::vector<gsMatrix<T> > > param_loc;

    bool exp_sampling = true;

    gsFittingParam<T> fitting_param_face(fitting_param);
    fitting_param_face.interiorKnots =
        fitting_param.interiorKnots_boundary;

    int ind = 0;
    for(unsigned i = 0;i < s;i++)
    {
        /*   if(i == 1 || i == 3)
             {*/
        fitting_param_face.output = fitting_param.output + "_face_"
            + util::to_string(i);

        gsTrimSurface<T> *trim_temp = faces_temp[i]->surf;
        gsTrimSurface<T> *trim_geom = faces_geom[i]->surf;

        gsTensorBSplineBasis<2> basis_loc
            = getBasisFromTemplate<2, T>
            (*trim_temp->getTP(),
             fitting_param.interiorKnots_boundary,
             fitting_param.degree);
        _result = fitFaceFixedParam
            (trim_temp, trim_geom, basis_loc,
             nPointsSample, fitting_param_face);

        GISMO_ASSERT(_result != NULL,
                     "Error during the parametrization of one face");
        exportGeometry<gsGeometry<T> >(_result,
                                       fitting_param_face.output);
        result = static_cast<gsGeometry<T>*>(_result);

        points_loc.push_back(std::vector< gsMatrix<T> >(1));
        param_loc.push_back(std::vector< gsMatrix<T> >(1));
        sample_boundary_points<T>(basis_loc, *trim_geom->getTP(),
                                  *trim_temp->getTP(),
                                  points_loc[ind], param_loc[ind], result);

        if(result != NULL)
            delete result;
        ind++;
        //   }
    }
    gsInfo << "Parametrization of the boundary finished"
           << std::endl;

    std::vector< gsMatrix<T> > points;
    std::vector< gsMatrix<T> > param;
    gsMatrix<T> bbox(3, 2);
    junctionPoints<T>(points_loc, param_loc, points, param);
    getBoundingBoxPts<T>(param[0], bbox);

    if(exp_sampling)
    {
        std::string path = fitting_param.output + "_pts";
        std::string path2 = fitting_param.output + "_image";
        gsWriteParaviewPoints(points[0], path2);
        gsWriteParaviewPoints(param[0], path);
    }
    gsTensorBSplineBasis<3> basis
        = getBasisFromBoundingBox<3, T>
        (bbox, fitting_param.interiorKnots,
         fitting_param.degree);
    fitting_param.deform_min = true;
    fitting_param.coeff_NL_global = 1.;
    fitting_param.coeff_linear_global = 1.;
    fitting_param.export_iterations = true;
    fitting_param.export_initialization = true;
    fitting_param.print_messages = true;
    fitting_param.num_lines = 5;
    fitting_param.coeff_tikhonov = 3000;
    fitting_param.inf_smoothing_parameter = 0.000001;

    gsInfo << "3D FITTING" << std::endl;
    return fitPoints<3, gsBasis<T>, 2, T>
        (points, param, basis, fitting_param);
}


template<class T> gsFunctionSet<T>*
fitFaceFixedParam(gsTrimSurface<T>* trim_temp,
                  gsTrimSurface<T>* trim_geom,
                  gsBasis<T>& basis,
                  int nPointsSample,
                  gsFittingParam<T> &fitting_param)
{
    std::vector< gsMatrix<T> > pts_loc;
    std::vector< gsMatrix<T> > param_loc;

    std::vector< gsMatrix<T> > pts, param;

    unsigned s_temp;

    s_temp = trim_temp->domain().numLoops();

    /// Sample the boundary
    for(unsigned j = 0;j < s_temp;j++)
    {
        pts_loc.push_back(gsMatrix<T>());
        param_loc.push_back(gsMatrix<T>());
        //  trim_geom->sampleLoop_into(j, nPointsSample, pts_loc[j]);
        trim_geom->domain().sampleLoop_into
            (j, nPointsSample, 2, pts_loc[j]);
        trim_temp->domain().sampleLoop_into
            (j, nPointsSample, 2, param_loc[j]);
    }

    GISMO_ENSURE(s_temp > 0, "ERROR: no trimming curve in this face");
    junctionPoints<T>(pts_loc, param_loc, pts, param);

    return fitPoints<2, gsBasis<T>, 2, T>
        (pts, param, basis, fitting_param);
}


template<short_t d, class T> gsTensorBSplineBasis<d>
getBasisFromTemplate(gsGeometry<T>& surf,
                     int nbInteriorKnots, int degree)
{
    /// compute the bounding box for the domain
    gsMatrix<T> bbox = surf.basis().support();
    return getBasisFromBoundingBox<d, T>(bbox, nbInteriorKnots,
                                   degree);
}

template<short_t d, class T> gsTensorBSplineBasis<d>
getBasisFromBoundingBox(gsMatrix<T> bbox,
                        int nbInteriorKnots, int degree)
{
    std::vector< gsKnotVector<T> > KV;
    for(unsigned i = 0; i < d; i++)
    {
        KV.push_back(gsKnotVector<T>
                     ( bbox(i, 0), bbox(i, 1),
                       nbInteriorKnots, degree + 1 ));
    }
    return gsTensorBSplineBasis<d>( KV );
}

template<short_t d, class T> gsTensorBSplineBasis<d>
getBasisFromBoundingBox(gsBoundingBox<T>& bbox,
                        int nbInteriorKnots, int degree)
{
    std::vector< gsKnotVector<T> > KV;
    for(unsigned i = 0; i < d; i++)
    {
        KV.push_back(gsKnotVector<T>
                     ( bbox.pMin(i, 0), bbox.pMax(i, 0),
                       nbInteriorKnots, degree + 1 ));
        gsInfo << util::to_string(bbox.pMin(i, 0)) << ";"
               << util::to_string(bbox.pMax(i, 0)) << "   ";

    }
    gsInfo <<  std::endl;
    return gsTensorBSplineBasis<d>( KV );
}


template<short_t d, class T> void
set_strong_bc_face(gsTrimSurface<T>& surf,
                   gsMultiBasis<T>& basis_face,
                   gsDofMapper mapper, int ind)
{
    GISMO_ASSERT(surf->domain().numLoops() == 1,
                 "The case where the surface has holes is not considered");
    std::vector< gsCurve<T> *>& loop = surf.loop(0).curves();
    unsigned nLoops = loop.size;
    GISMO_ASSERT(basis_face.nBases() == nLoops,
                 "We only consider the mid-edge splitting: the number of patches should be the same as the number of edges for each patch");

    for(unsigned i = 0;i < nLoops;i++)
    {
        /// We project the boundary on the basis
        int iM = (i + nLoops - 1) % nLoops;
        gsCurve<T> *cur_curv = loop[i];
        gsCurve<T> *curvM = loop[iM];
        gsBasis<T>& _tensor_b = basis_face.basis(i);
        gsTensorBasis<d, T>& tensor_b = static_cast<gsTensorBasis<d, T>&>(_tensor_b);

        /// y coordinates
        gsMatrix<index_t> indices_y = tensor_b.boundary(1);
        /// x coordinates
        gsMatrix<index_t> indices_x = tensor_b.boundary(0);
        gsBasis<T>& b1 = tensor_b.component(0);
        gsBasis<T>& b2 = tensor_b.component(1);
    }

}

/// The curve is split into two components
/// begin is true if the first part of the curve is considered
/// (the second part is considered otherwise)
template<class T> gsGeometry<T>*
get_proj_boundary(gsBasis<T>& new_basis, gsCurve<T>& curv,
                  bool begin)
{
    gsMatrix<T> supp_c = curv.basis()->support();
    if(begin)
        supp_c(0, 1) = supp_c(0, 0)
            + (supp_c(0, 1) - supp_c(0, 0)) * 0.5;
    else
        supp_c(0, 0) = supp_c(0, 1)
            - (supp_c(0, 1) - supp_c(0, 0)) * 0.5;
    gsMatrix<T> supp_b = new_basis->support();

    short_t targetDim = curv.targetDim();

    gsMatrix<T> anch = new_basis.anchors();
    int nAnch = anch.rows();
    gsMatrix<T> image(targetDim, nAnch);
    gsMatrix<T> param(targetDim, nAnch);

    T diff = supp_c(0,0) - supp_b(0,0);
    T coeff = ( supp_c(0,1) - supp_c(0,0) )
        / ( supp_b(0,1) - supp_b(0,0) );
    for(unsigned i = 0;i < nAnch;i++)
        param(i,0) = diff + coeff * anch(i,0);
    curv.eval_into(param, image);
    return interpolateData(image, anch).release();
}

} /// namespace gismo

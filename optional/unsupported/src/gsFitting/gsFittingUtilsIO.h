/** @file gsFittingIOUtils.h

    @brief This file contains various input and output routines used for fitting

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsGeometrySlice.h>

namespace gismo
{

/// Exports the result of a fitting
template<class GenGeom, class T>
void exportGeometry(gsFunctionSet<T>* map,
                    std::string output,
                    bool exp_xml = false,
                    int NparameterLines = 50,
                    int dir_slices = 2, int nSlices = 5);

/// Exports some slices of the geometry (only for dim=3)
template<class T>
void exportSlices(gsGeometry<T>& map, std::string output,
                  int dir_slices, int nSlices, int nPoints);


/// Returns the number of patches
template<class T>
int getPoints(std::vector< gsMatrix<T> >& params,
              std::vector< gsMatrix<T> >& points,
              const std::string& filename);


/* Read the two files given by the user. The input can be:
   - two single or multipatch geometries being the template and the target geometries.
   The image of the template domain must, in that case, be the target geometry.
   - two set of points (p_i) and (q_i):
   the image of q_i must be p_i in that case.
   - two gsSolid being the parameter domain and the surface to parameterize.
   We first parametrize each face of the solid in that case. Then, we construct the 3D
   parametrization of the volume such that the image of a face of the parameter domain
   is it's associated face in the surface.
   In the case the topology of the multipatch must be given, it is contained in fic_topo.
*/
template<class T>
int readInput(gsFittingParam<T> &fitting_param);

/// taken from motor/jku/gsMotorIOUtils
void writeKnots(const gsGeometry<>& map,
                const std::string& outputPrefix);
void writeKnots(const gsGeometry<>& map,
                const int id,
                const std::string& outputPrefix);

} /// namespace gismo

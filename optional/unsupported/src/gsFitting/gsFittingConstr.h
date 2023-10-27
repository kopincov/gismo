/** @file gsFittingConstr.h

    @brief Files containing all the routines permitting
    to construct the different objects necessary to use gsFitting.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingIterative.h>
#include <gsFitting/gsFittingHierar.h>
#include <gsFitting/gsFittingReinit.h>
#include <gsFitting/gsFittingIntegrandNL.h>

#include <gsFitting/gsFittingEnergy.h>
#include <gsFitting/gsFittingBase.h>
#include <gsFitting/gsLeastSquaresLMult.h>
#include <gsFitting/gsFittingEnergy.h>

#include <gsFitting/gsBoundariesCoupling.h>
#include <gsFitting/gsFittingGeometry.h>

namespace gismo
{

/**
   Interface containing a function that initializes
   a gsFitting. This permits to use only one fitting routine
   (called fitData) for any given data to be fitted
 */
template<short_t d, class GenBasis, class T>
class gsFittingInitializer;


/**
   Class containing a function that initializes
   a gsFitting in the case where the least squares method is used.
 */
template<short_t d, class GenBasis, class T>
class gsFittingInitializerLS;

/**
   Class containing a function that initializes
   a gsFitting in the case where the L2 distance minimization is used.
 */
template<short_t d, class GenBasisDom, class GenGeomBC,
         unsigned o, class T>
class gsFittingInitializerContinuous;

/// Remove the geometries
template<short_t d, class GenGeom, unsigned o, class T>
void remove_strongBC(GenGeom& temp, GenGeom& geom,
                     std::vector<int> ind);


/// Creates the quadrature associated with the smoothing
/// (adds all the linear and nonlinear energies)
template<short_t d, class GenBasis,
         unsigned o, class T> void
addSmoothingQuadrature(GenBasis& basis, int dim_im,
                       gsFittingQuadrEnergy<GenBasis, T>& ener,
                       gsFittingParam<T>& param);

template<short_t d, class GenBasis, unsigned o, class T> void
addSmoothingIntegrands(gsFittingParam<T>& param,
                       gsFittingQuadrature<d, GenBasis, o, T>&
                       quadr);

/// Adds the nonlinear smoothing energy to the vector of integrands integr
template<short_t d, class GenBasis, unsigned o, class T> void
add_NL_smoothing(gsFittingQuadrature<d, GenBasis, o, T>& quadr,
                 gsFittingParam<T>& param);


/// Creates the quadrature associated with the geometries to be fitted
template<short_t d, class GenBasisDom,
         class GenGeomBC, unsigned o, class T> void
addBCQuadrature(GenBasisDom& basis, GenGeomBC& domain,
                GenGeomBC& geom, unsigned dim_im,
                gsFittingQuadrEnergy<GenBasisDom, T>& ener,
                int size_gauss_quadr = 0);

/// Adds the winslow energy to the vector of integrands integr
/// The cases 2d and 3d must be distinguished.
template<class GenBasis, class T> void
add_winslow_smoothing(gsFittingQuadrature<2, GenBasis, 2, T>&
                      quadr, T coeff_NL, T coeff_winslow);
template<class GenBasis, class T> void
add_winslow_smoothing(gsFittingQuadrature<3, GenBasis, 2, T>&
                      quadr, T coeff_NL, T coeff_winslow);

/// Constructs the boundary conditions associated
/// with the border and the given indices.
template<class T> gsBoundaryConditions<T>
set_fixed_border(gsMultiPatch<T>& border, std::vector<int> ind);

/// Returns the indices of the geometries that defines
/// Dirichlet boundary conditions
std::vector<int> indices_bc(int dim, const std::vector<int>& indices);

/** Constructs a mapping minimizing a smoothing energy
 and the distance to some given points representing
 the template and the target geometries.
 These geometries should be given in case they are explicitly known.
 This permits to resample in case we refine the basis.
 This routine calls fitData.
*/
template<short_t d, class GenBasis,
         unsigned o, class T> gsFunctionSet<T>*
fitPoints(std::vector< gsMatrix<T> >& points,
          std::vector< gsMatrix<T> >& parameters,
          GenBasis& basis, gsFittingParam<T>& param,
          gsFunctionSet<T>* templ = NULL,
          gsFunctionSet<T>* geom = NULL);


/**
   Constructs a mapping minimizing a smoothing energy
 and the distance to a given target geometry, i.e.,
 the distance between the image of the geometry temp and
 the geometry geom.
 This routine calls fitData.
*/
template<short_t d, class GenBasisDom, class GenGeomBC,
         unsigned o, class T> gsFunctionSet<T>*
fitGeometry(GenGeomBC* geom, GenGeomBC* temp,
            GenBasisDom& basisDom, gsFittingParam<T>& param);

/**
   Generic routine called to fit some data.
 */
template<short_t d, class GenBasis, class T> gsFunctionSet<T>*
fitData(int dim_im, GenBasis& basis, gsFittingParam<T>& param,
        gsFittingInitializer<d, GenBasis, T>& initializer);

template<short_t d, class GenBasis, class T> gsFunctionSet<T>*
fitData(int dim_im, GenBasis& basis, gsFittingParam<T>& param,
        gsFittingGeometryCoupling<T>& coupling);

/**
   Interface containing a function that initializes
   a gsFitting. This permits to use only one fitting routine
   (called fitData) for any given data to be fitted
 */
template<short_t d, class GenBasisDom, class T>
class gsFittingInitializer
{
public:
    gsFittingInitializer()  {  }
    virtual void init(gsFittingIter<d, GenBasisDom, T>& fitting)
    {  GISMO_NO_IMPLEMENTATION;  }
};


/**
   Class containing a function that initializes
   a gsFitting in the case where the least squares method is used.
 */
template<short_t d, class GenBasis, class T>
class gsFittingInitializerLS
    : public gsFittingInitializer<d, GenBasis, T>
{
    typedef gsFittingInitializer<d, GenBasis, T> Base;
public:
    /** TODO: set the template and the geometry
        to resample the points if necessary
    */
    gsFittingInitializerLS(std::vector< gsMatrix<T> >& points,
                           std::vector< gsMatrix<T> >& parameters,
                           gsFunctionSet<T>* templ,
                           gsFunctionSet<T>* geom)
    : Base::gsFittingInitializer(), m_points(points),
    m_parameters(parameters)
    {
        m_templ = templ;
        m_geom = geom;
    }

    void init(gsFittingIter<d, GenBasis, T>& fitting)
    {
        fitting.energy().
            resetPointsLS(m_points, m_parameters);
    }

private:
    std::vector< gsMatrix<T> >& m_points;
    std::vector< gsMatrix<T> >& m_parameters;
    gsFunctionSet<T>* m_templ;
    gsFunctionSet<T>* m_geom;
};


/**
   Class containing a function that initializes
   a gsFitting in the case where the L2 distance minimization is used.
 */
template<short_t d, class GenBasisDom, class GenGeomBC,
         unsigned o, class T>
class gsFittingInitializerContinuous
    : public gsFittingInitializer<d, GenBasisDom, T>
{
    typedef gsFittingInitializer<d, GenBasisDom, T> Base;
public:
    gsFittingInitializerContinuous(GenGeomBC* geom,
                                   GenGeomBC* temp)
    : Base::gsFittingInitializer()
    {
        m_geom = geom;
        m_temp = temp;
    }

    ~gsFittingInitializerContinuous()
    {
    }

    void init(gsFittingIter<d, GenBasisDom, T>& fitting)
    {
        GISMO_ASSERT(m_temp != NULL && m_geom != NULL,
                     "Should not be NULL");

        refineBoundary(*m_temp, fitting.basis());
        addBCQuadrature<d, GenBasisDom, GenGeomBC, o, T>
            (fitting.basis(), m_temp, m_geom, fitting.dim_im(),
             fitting.quadrEner());
    }

private:
    GenGeomBC* m_geom;
    GenGeomBC* m_temp;
};


} /// namespace gismo

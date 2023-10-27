/** @file gsBoundariesCoupling.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsFitting/gsFittingLinSyst.h>
#include <gsFitting/gsFittingQuadrature.h>
#include <gsFitting/gsFittingIntegrandNL.h>

namespace gismo
{

/// @brief In case a border of the template and target geometries
/// are defined by trimmed surfaces.
/// A mapping must be constructed between the trimmed domains of
/// these two trimmed surface. This class permits to construct
/// this mapping
template<class T>
class gsMapTrimDom
{

    /// In case we fit the
    gsFunctionSet<T>* m_param_trimmed;

    ~gsMapTrimDom()
    {
        if(m_param_trimmed != NULL)
            delete m_param_trimmed;
    }
};

struct type_border{
public:
    enum type
    {
        unknownType = -1,
        pts         =  0, ///< Defined by points
        geom        =  1, ///< Defined by geometries
        trim_surf   =  2, ///< Defined by trimming surfaces
        side        =  3  ///< Defined by a side of a multi-tensor domain
    };
};

struct type_coupling{
public:
    enum type
    {
        unknownType = -1,
        PtsPts      =  0, ///< Defined by two sets ofpoints
        GeomGeom    =  1, ///< Defined by two geometries
        TrimTrim    =  2, ///< Defined by two trimming surfaces
        StrongPts   =  3,  ///< Strong bc given by points
        StrongGeom  =  4  ///< Strong bc given by a geometry
    };
};

struct SideFitting{
    bool sgn_orien;
    patchSide side;

    SideFitting(patchSide _side, bool _sgn_orien)
    {
        sgn_orien = _sgn_orien;
        side = _side;
    }
    SideFitting(){  }
};


/* @brief Class containing the different types
   of data defining a border used for fitting.
   These Borders are necessarily included in one patch
*/
template<class T>
class gsFittingBorder
{
protected:

    /// The delimiting points
    gsMatrix<T> m_points;
    /// The delimiting geometry
    gsGeometry<T>* m_geometry;

    /// The trimmed surface
    gsTrimSurface<T> m_trimmed;

    SideFitting m_strong;

    /// The type of the border
    type_border::type m_type;


public:
    gsFittingBorder()
    : m_type(type_border::unknownType) {  init();  }

    gsFittingBorder(gsMatrix<T>& points)
    : m_points(points), m_type(type_border::pts) {  init();  }

    gsFittingBorder(gsGeometry<T>& geom)
    : m_type(type_border::geom)
    {
        init();
        m_geometry = &geom;
    }

    gsFittingBorder(gsTrimSurface<T>& trimmed)
    : m_trimmed(trimmed), m_type(type_border::trim_surf) {  init();  }

    /*gsFittingBorder(index_t patch, boundary::side s)
      : m_strong(patch, s), m_type(type_border::side) {  init();  }*/

    gsFittingBorder(const gsFittingBorder<T>& source);


    gsMatrix<T>& points() {  return m_points;  }
    gsGeometry<T>& geometry()
    {
        GISMO_ASSERT(m_geometry != NULL,
                     "The geometry is NULL");
        return *m_geometry;
    }
    gsTrimSurface<T>& trimmed() {   return m_trimmed;  }
    SideFitting& strong() {   return m_strong;  }
    type_border::type& type(){  return m_type;  }

    void getBoundingBox(gsMatrix<T>& bbox);
    index_t dim_im();

private:
    void init(){  m_geometry = NULL;   }
};



template<class T>
class gsFittingBorderCoupling
{
protected:

    gsFittingBorder<T>& m_cond_template;
    gsFittingBorder<T>& m_cond_target;

    type_coupling::type m_type;

public:

    gsFittingBorderCoupling(gsFittingBorder<T>& cond_template,
                            gsFittingBorder<T>& cond_target);

    gsFittingBorderCoupling(const gsFittingBorderCoupling<T>& source)
    : m_cond_template(source.m_cond_template),
      m_cond_target(source.m_cond_target), m_type(source.m_type) {  }

    template<class GenBasis>
    void associatePts(gsFittingEnergy<GenBasis, T>& energy,
                      GenBasis& basis, index_t ind);

    template<class GenBasis>
    void associateTrim(gsFittingEnergy<GenBasis, T>& energy,
                       GenBasis& basis);

    template<short_t d, class GenBasis>
    void associateGeom(gsFittingEnergy<GenBasis, T>& energy,
                       GenBasis& basis, index_t ind);

    template<short_t d, class GenBasis>
    void associate(gsFittingBase<d, GenBasis, T>& fitting,
                   index_t ind);

    template<short_t d, class GenBasis>
    void associateStrong(gsFittingBase<d, GenBasis, T>& fitting);

    /// Construct the data used for the association between
    /// the two borders (for example construct the parametrization
    /// of the trimming curve if needed, sample points from
    /// the curve if needed, ...)
    template<class GenBasis> void
    initialize(gsFittingParam<T>& param, GenBasis& basis);

    /// Construct a parametrization of the trimmed surface
    /// in order to add it to the gsFitting
    template<class GenBasis> void
    initializeTrim(gsFittingParam<T>& param, GenBasis& basis);

#if ! EIGEN_HAS_RVALUE_REFERENCES
    /// \brief Swap with another gsFittingBorderCoupling
    void swap(gsFittingBorderCoupling& other)
    {
        std::swap(m_cond_target, other.m_cond_target);
        std::swap(m_cond_template, other.m_cond_template);
        std::swap(m_type, other.m_type);
    }

    /// Assignment operator (uses copy-and-swap idiom)
    gsFittingBorderCoupling& operator= ( gsFittingBorderCoupling other )
    {
        this->swap(other);
        return *this;
    }
#endif

}; /// gsFittingBorderCoupling

} /// namespace gismo

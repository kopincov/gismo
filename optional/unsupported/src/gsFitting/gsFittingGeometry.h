/** @file gsFittingGeometry.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

namespace gismo
{

/**
   @brief Definition of a geometry in gsFitting.
   A geometry can either be:
   - a multi-tensor domain, in which case each tensor domain is defined by a bounding box
   - a manifold that is defined using boundary representations.
   In the context of gsFitting, the border can be given by a set of points, a geometry or a trimmed geometry
   Note that the trimmed geometry will first be parametrized before
   calling the fitting
*/
template<class T>
class gsFittingGeometry
{
protected:
    typedef std::vector< gsFittingBorder<T> > ContainerPatch;

    /// The boundaries of the geometry
    std::vector<ContainerPatch> m_boundaries;
    std::vector< gsMatrix<T> > m_bbox;

public:
    gsFittingGeometry(index_t nPatches = 0)
    : m_boundaries(nPatches) { }

    gsFittingGeometry(std::vector<ContainerPatch>& boundaries)
    : m_boundaries(boundaries) { }

    gsFittingGeometry(std::vector< gsMatrix<T> >& pts,
                      bool unique_patch);
    gsFittingGeometry(gsMultiPatch<T>& geom, bool unique_patch);

    ////////// The geometry has only one patch ////////////
    gsFittingGeometry(ContainerPatch& boundaries)
    : m_boundaries(boundaries) { }

    gsFittingGeometry(gsMatrix<T>& pts)
    : m_boundaries()
    {
        m_boundaries.push_back
            ( std::vector< gsFittingBorder<T> >
              (gsFittingBorder<T> (pts) ) );
    }

    gsFittingGeometry(gsGeometry<T>& geom)
    : m_boundaries()
    {
        m_boundaries.push_back
            ( std::vector< gsFittingBorder<T> >( ) );
        m_boundaries[0].push_back
            (gsFittingBorder<T> (geom));
    }

    /// Constructs a geometry from a solid
    gsFittingGeometry(gsSolid<T>& solid);

    /// Construct a geometry from a trimmed domain in the plane.
    /// Note that here the trimmed domain is not the border.
    /// The border is given by the trimming curves.
    gsFittingGeometry(gsPlanarDomain<T>& domain);

    gsFittingGeometry(const gsFittingGeometry<T>& source)
    : m_boundaries(source.m_boundaries),
      m_bbox(source.m_bbox) {  }

    gsFittingGeometry<T>& operator=(const gsFittingGeometry<T>& source)
    {
        m_bbox = source.m_bbox;
        m_boundaries = source.m_boundaries;
        return *this;
    }

    //////////////////////////////////////////////////////

    index_t nPatches()
    {
        if(m_boundaries.size() > 0 )
            return m_boundaries.size();
        else
            return m_bbox.size();
    }

    void add_boundary(index_t patch, gsFittingBorder<T>& boundary)
    {
        GISMO_ASSERT(patch < m_boundaries.size(), "Number of patch lower that the index given");
        m_boundaries[patch].push_back(boundary);
    }

    index_t nBorder(index_t i)
    {
        GISMO_ASSERT(i < nPatches(), "Error in the index");
        if(m_boundaries.size() > 0)
            return m_boundaries[i].size();
        return 0;
    }

    index_t dim_im()
    {
        GISMO_ASSERT(m_boundaries.size() > 0, "No border!!!");
        GISMO_ASSERT(m_boundaries[0].size() > 0,
                     "The patch has no border!!!");
        return m_boundaries[0][0].dim_im();
    }

    void getBoundingBox(gsMatrix<T>& res);

    gsFittingBorder<T>& border(index_t patch, index_t pos)
    {
        GISMO_ASSERT(patch < (index_t)m_boundaries.size()
                     && pos < (index_t)m_boundaries[patch].size(),
                     "Indices invalid");
        return m_boundaries[patch][pos];
    }
};

template<class T>
class gsFittingGeometryCoupling
{
protected:
    typedef std::vector< gsFittingBorderCoupling<T> > ContainerPatch;

    /// The template geometry
    gsFittingGeometry<T>* m_template;

    /// The target geometry
    gsFittingGeometry<T>* m_target;

    /// All the coupling between the template and the target
    std::vector<ContainerPatch> m_coupling;

public:
    gsFittingGeometryCoupling(gsFittingGeometry<T>& templ,
                              gsFittingGeometry<T>& target,
                              std::vector<ContainerPatch>&
                              coupling)
    : m_template(&templ), m_target(&target),
      m_coupling(coupling) { }

    gsFittingGeometryCoupling(const gsFittingGeometryCoupling
                              <T>& source)
    : m_template(source.m_template), m_target(source.m_target),
      m_coupling(source.m_coupling) { }

    /// We associate (in the same order) all the boundaries
    /// of the template with the boundaries of the target
    gsFittingGeometryCoupling(gsFittingGeometry<T>& templ,
                              gsFittingGeometry<T>& targ);

    gsFittingGeometryCoupling<T>&
    operator=(const gsFittingGeometryCoupling<T>& source)
    {
        m_template = source.m_template;
        m_target = source.m_target;
        m_coupling.clear();
        unsigned s = source.m_coupling.size();
        for(unsigned i = 0;i < s;i++)
            m_coupling.push_back(source.m_coupling[i]);
        return *this;
    }

    /// Compute a bounding box for each patch of the template.
    /// Then, constructs a basis on this bounding box.
    template<unsigned d> gsTensorBSplineBasis<d, T>
    getBasisTrimming(gsFittingParam<T>& param_fitting);

    /// Split the boundary conditions and
    void split_bc(std::vector<index_t> indices_split);

    ///
    template<class GenBasis> void
    initialize(gsFittingParam<T>& param, GenBasis& basis)
    {
        unsigned s1 = m_coupling.size();
        for(unsigned i = 0;i < s1;i++)
        {
            unsigned s2 = m_coupling[i].size();
            for(unsigned j = 0;j < s2;j++)
                m_coupling[i][j].template
                    initialize<GenBasis>(param, basis);
        }
    }

    template<short_t d, class GenBasis>
    void associate(gsFittingBase<d, GenBasis, T>& fitting)
    {
        unsigned s = m_coupling.size();
        for(unsigned i = 0;i < s;i++)
        {
            unsigned s2 = m_coupling[i].size();
            for(unsigned j = 0;j < s2;j++)
                m_coupling[i][j].template
                    associate<d, GenBasis>(fitting, i);
        }
    }

};

} /// namespace gismo

/** @file gsSelector.h

    @brief Chooses which lower level basis to call in a given point.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
**/

#pragma once

#include <fstream>

#include <gsCore/gsExport.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunctionSet.h>

#include <gsRemappedBasis/tree.h>
#include <gsRemappedBasis/gsBoxList.h>
#include <gsRemappedBasis/gsDomainMap.h>

namespace gismo
{

/** A data structure that allows to determine in which subdomain a
 * given point is. It is used in gsRemappedBasis and described in
 * \cite bm2016, Section 2.1, where it is called partition
 * \f$\mathbf{D}\f$.
 *
 * Each patch (gsSelector can represent multipatch domains) is
 * represented by a gsDomainMap.
 *
 * The gsSelector can be constructed from a gsBoxList. Cf. the
 * detailed description in gsDomainMap.*/
class GISMO_EXPORT gsSelector
{
    typedef unsigned basisIdT;
    typedef int directionT;
    typedef size_t  patchIdT;
    
public:
    typedef gsDomainMap::NodeId NodeId;

public: // Initializers

    void initFromBoxes (
        const gsBoxList        &boxes,
        const gsMatrix<real_t> &boundingBox,
        size_t patch = 0);

    void initFromBoxesMax (
        const gsBoxList        &boxes,
        const gsMatrix<real_t> &boundingBox,
        size_t patch = 0);

    void initFromBoxesMin (
        const gsBoxList        &boxes,
        const gsMatrix<real_t> &boundingBox,
        size_t patch = 0);

    void addBoxesMin(
        const gsBoxList        &boxes,
        size_t patch = 0);

    void addBoxesMax(
        const gsBoxList        &boxes,
        size_t patch = 0);

public: // Getters

    gsDomainMap&       patch(size_t p);
    const gsDomainMap& patch(size_t p) const;

    size_t             patchNum () const;

    /**
     * Returns the basisId of the basis capable of depicting all the basis functions
     * active in the point @a p.
     */
    template <typename PointType>
    basisIdT getBasisAt (const PointType &p, patchIdT patch=0) const
    {
        return m_domainMap[patch].getBasisAt(p);
    }

    /**
     * Returns the basisId of the basis capable of depicting all the basis functions
     * active in the given points.
     * @param points , each given as a column of the matrix;
     * @param out returned basisId's (one per point).
     */
    void getBasisAt ( const gsMatrix<real_t> &points,
                      std::vector<basisIdT>  &out,
                      patchIdT               patch=0) const;

public: // Exports

    gsBoxList asBoxList( size_t patch=0 ) const;

    void exportToTex(std::string filename) const;

private: // members
    std::vector<gsDomainMap> m_domainMap;
};

} // namespace gismo

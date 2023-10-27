/** @file gsIteratorFitting.h

    @brief Contains iterators used for gsFitting:
     gsPointIterator: iterator on gsPointContainer
     gsGenDomainIterator: iterator on domains that can be either
   single patch or multipatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Y. Masson
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsFitting/gsPointContainer.h>

namespace gismo
{

/* @brief
   Class gsPointIterator: iterator on gsPointContainer. These points can be contained in a single patch or a multipatch geometry.
*/
template<class T>
class gsPointIterator
{

public:
    /// constructor
    gsPointIterator(gsPointContainer<T>& container) :
    m_container(container)
    {
        reset();
    }

    /// Destructor
    ~gsPointIterator(){  }

    /// Returns the current patch
    inline index_t patch(){ return m_patch; }

    /// Returns the current position in the current patch
    inline index_t pos(){ return m_pos; }

    /// Goes to the next position. Return true if the iterator is still active
    /// (not the end of the iteration)
    inline bool next()
    {
        if(m_is_good)
        {
            m_pos++;
            if(m_pos == m_container.size(m_patch))
            {
                m_pos = 0;
                m_patch++;
            }
        }
        isGoodActualization();
        return m_is_good;
    }


    /// Returns true if the iterator is still active
    /// (not the end of the iteration)
    inline bool good(){ return m_is_good; }

    /// Returns the point at the current position if m_is_good.
    inline typename gsMatrix<T>::Column currPoint ()
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_container.point(m_patch, m_pos);
    }

    /// Returns the parameter at the current position if m_is_good
    inline typename gsMatrix<T>::Column currParam ()
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_container.param(m_patch, m_pos);
    }

    /// Returns the basis at the current patch if m_is_good
    inline gsBasis<T>& currentBasis(gsMultiBasis<T> &basis)
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return basis[m_patch];
    }

    /// Returns the basis at the current patch if m_is_good
    inline gsBasis<T>& currentBasis(gsBasis<T> &basis)
    { return basis;  }

    /// Returns the global index of the current position if m_is_good
    inline int globalIndex()
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_container.globalIndex(m_patch, m_pos);
    }

    /// Returns the point of the set pts at the position of the iterator
    /// if m_is_good
    inline typename gsMatrix<T>::Column
    getPointUsingIterator(std::vector< gsMatrix<T> >& pts)
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return pts[m_patch].col(m_pos);
    }

    /// Sets some point in the set pts at the position of the iterator
    /// if m_is_good
    inline void setPointUsingIterator(std::vector< gsMatrix<T> >& pts,
                                      typename gsMatrix<T>::Column& pt)
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        pts[m_patch].col(m_pos) = pt;
    }

private:
    /// Sets the iterator to the initial position
    bool reset()
    {
        m_pos = 0;
        m_patch = 0;
        isGoodActualization();
        return m_is_good;
    }

    /// Actualize the value of m_is_good (false if the iteration is over)
    inline void isGoodActualization(){
        m_is_good = m_patch < m_container.nPatches()
                              && m_pos < m_container.size(m_patch);
    }

    /// The point container on which is defined the iterator
    gsPointContainer<T>& m_container;

    /// the position of the point iterator (point set)
    index_t m_pos;
    /// the current patch
    index_t m_patch;
    /// is the iteration over
    bool m_is_good;

}; // class gsPointIterator


/**
   @brief
   Class gsGenDomainIterator: an iterator on domains that can be either
   single patch or multipatch.

**/
template<class T>
class gsGenDomainIterator
{


public:
    /// Constructors

    /// Creates a iterator on the domain of a gsBasis
    gsGenDomainIterator(gsBasis<T> &basis){ update(basis); }

    /// Creates a iterator on the domain of a gsMultiBasis
    gsGenDomainIterator
    (gsMultiBasis<T> &basis){ update(basis); }

    /// Creates a iterator on the domain of a gsMultiPatch
    gsGenDomainIterator
    (gsMultiPatch<T> &patches){ update(patches); }

    /// Copies a gsGenDomainIterator
    gsGenDomainIterator(gsGenDomainIterator<T>& it)
    : m_patch(it.m_patch), m_dom_iter(it.m_dom_iter)
    {
        m_total_numElements = it.m_total_numElements;
        reset();
    }

    /// Creates an empty iterator
    gsGenDomainIterator(){ }

    /// Destructor
    ~gsGenDomainIterator(){ }

    /// Resets the iterator to initial position on the domain
    /// of a gsBasis
    void setBasis(gsBasis<T> &basis){ update(basis); }

    /// Resets the iterator to initial position on the domain
    /// of a gsMultiBasis
    void setBasis(gsMultiBasis<T> &basis){ update(basis); }

    /// Resets the iterator to initial position on the domain
    /// of a gsMultiPatch
    void setBasis(gsMultiPatch<T> &patches){ update(patches); }

    /// Sets (resets) the number of elementary cells
    /// in the domain
    void init_numElements()
    {
        unsigned s = m_dom_iter.size();
        m_total_numElements = 0;
        for(unsigned i = 0;i<s;i++)
        {
            m_numElements.push_back
                (m_dom_iter[i]->numElements());
            m_total_numElements += m_numElements[i];

        }
    }

    /// Returns the dimension of the parameter domain
    inline int paramDim()
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_dom_iter[0]->dim();
    }

    /// Goes to the next elementary cell. Returns true if
    /// the iteration is not over
    inline bool next()
    {
        if(m_is_good)
        {
            m_dom_iter[m_patch]->next();
            if(! m_dom_iter[m_patch]->good())
            {
                m_patch++;
                if((unsigned) m_patch >= m_dom_iter.size())
                {
                    m_is_good = 0;
                    return false;
                }
                else
                {
                    isGoodActualization();
                    return m_is_good;
                }
            }
            return true;
        }
        return false;
    }

    /// Actualize the value of m_is_good (false if the iteration is over)
    inline void isGoodActualization()
    {
        if((unsigned) m_patch >= m_dom_iter.size())
            m_is_good = false;
        else
            m_is_good = m_dom_iter[m_patch]->good();
    }

    /// Returns true if the iterator is still active
    /// (not the end of the iteration)
    inline bool good(){ return m_is_good; }

    /// Returns the lower corner of the current elementary cell
    inline const gsVector<T>& lowerCorner() const
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_dom_iter[m_patch]->lowerCorner();
    }

    /// Returns the upper corner of the current elementary cell
    inline const gsVector<T>& upperCorner() const
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_dom_iter[m_patch]->upperCorner();
    }

    /// Returns the central point of the current elementary cell
    inline const gsVector<T>& centerPoint() const
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_dom_iter[m_patch]->centerPoint();
    }

    /// Returns the basis of the current patch
    inline gsBasis<T>& currentBasis(gsBasis<T> &basis)
    {  return basis;   }

    /// Returns the basis of the current patch
    inline gsBasis<T>& currentBasis(gsMultiBasis<T> &basis)
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return basis[m_patch];
    }

    /// Returns the geometry associated with the current patch in geom
    inline gsGeometry<T>& currentGeom(gsMultiPatch<T> &geom)
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return geom.patch(m_patch);
    }

    /// Returns the geometry associated with the current patch in geom
    inline gsGeometry<T>& currentGeom(gsGeometry<T> &geom)
    { return geom;  }

    /// Returns the current patch
    inline int patch(){
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_patch;
    }

    /// Returns the total number of elementary cells
    inline int total_numElements()
    {   return m_total_numElements;   }

    /// Returns the number of elementary cells in some patch
    inline int numElement(int patch)
    {   return m_numElements[patch];   }

    /// Returns the number of patches
    inline int nPatches()
    { return m_numElements.size();  }


private:
    /// Sets the iterator to the initial position
    inline bool reset()
    {
        index_t size = m_dom_iter.size();
        m_patch = 0;
        if(size == 0){
            m_is_good = 0;
            return 0;
        }
        isGoodActualization();
        return m_is_good;
    }

    /// Initialize the domain iterator. This domain iterator
    /// will be defined on the domain of a gsBasis
    void update(gsBasis<T> &basis)
    {
        m_dom_iter.clear();
        m_dom_iter.push_back(basis.makeDomainIterator());
        init_numElements();
        reset();
    }

    /// Initialize the domain iterator. This domain iterator
    /// will be defined on the domain of a gsMultiBasis
    void update(gsMultiBasis<T> &basis)
    {
        m_dom_iter.clear();
        unsigned s = basis.nBases();
        for(unsigned i = 0;i<s;i++)
            m_dom_iter.push_back(basis[i].makeDomainIterator());
        init_numElements();
        reset();
    }

    /// Initialize the domain iterator. This domain iterator
    /// will be defined on the domain of a gsMultiPatch
    void update(gsMultiPatch<T> &patches)
    {
        m_dom_iter.clear();
        unsigned s = patches.nPatches();
        for(unsigned i = 0;i<s;i++)
            m_dom_iter.push_back(patches[i].basis()
                                 .makeDomainIterator());
        init_numElements();
        reset();
    }

protected:
    /// The current patch
    index_t m_patch;

    /// The total number of elementary cells
    index_t m_total_numElements;

    /// The number of elementary cells in each patch
    std::vector<index_t> m_numElements;

    /// The domain iterators on each patch
    std::vector<typename gsBasis<T>::domainIter> m_dom_iter;

    /// True if the iteration is not over
    bool m_is_good;


}; // class gsGenDomainIterator



/**
   @brief Class gsGenGeomIterator: an iterator on a geometry.
   This geometry can be either single patch or multipatch.
   This class is abstract. Three different cases are considered:
     - the geometry is the identity (this class). The functions called
       on the geometry are simply applied on the gsGenDomainIterator.
       Then, it can be single patch or multipatch depending on the domain
       of the iterator.
     - the geometry is single patch. The image of the geometry must then
       be a single patch in the current implementation.
       See class gsGenGeomIteratorSimplePatch.
     - the geometry is multipatch. The image of the geometry must then
       be a single patch in the current implementation.
       See class gsGenGeomIteratorMultiPatch.
   \ingroup Modeling
**/
template<class T>
class gsGenGeomIterator : public gsGenDomainIterator<T>
{

    using gsGenDomainIterator<T>::m_dom_iter;

public:

    /// constructors (see gsGenDomainIterator)
    gsGenGeomIterator(gsBasis<T>& basis) :
    gsGenDomainIterator<T>::gsGenDomainIterator(basis)
    {  }
    virtual ~gsGenGeomIterator(){  }
    gsGenGeomIterator(gsMultiBasis<T>& basis) :
    gsGenDomainIterator<T>::gsGenDomainIterator(basis)
    {  }
    gsGenGeomIterator(gsMultiPatch<T>& basis) :
    gsGenDomainIterator<T>::gsGenDomainIterator(basis)
    {  }

    /// Returns the upper corner of the image of the elementary cell
    virtual const gsVector<T> imageUpperCorner() const
    { return this->upperCorner(); }

    /// Returns the lower corner of the image of the elementary cell
    virtual const gsVector<T> imageLowerCorner() const
    { return this->lowerCorner(); }

    /// Returns the central point of the image of the elementary cell
    virtual const gsVector<T> imageCenterPoint() const
    { return this->centerPoint(); }

    /*
    /// res: cols() = number of points giving the directions (d+1)
    ///      rows() = dimension of the ambient space
    virtual void imageElementaryCell_into(gsMatrix<T> &res) const
    {
        unsigned d = m_dom_iter[0]->dim();
        res.col(0) = this->lowerCorner();
        gsVector<T> diag = this->upperCorner()
            - this->lowerCorner();

        for(unsigned i = 0;i < d;i++)
        {
            res.col(i + 1) = this->lowerCorner();
            res.col(i + 1).row(i) += diag.row(i);
        }
        }*/

    /// Evaluates the image of points on the current patch by the geometry
    virtual const void eval_into(gsMatrix<T>& points,
                                 gsMatrix<T>& res) const
    {  res = points;  }

    /// Returns the dimension of the domain on which is
    /// embedded the geometry
    virtual const unsigned basisDomDim() const
    {  return m_dom_iter[0]->dim();  }

    /// Returns the basis of the patch in which is included the image
    /// by the geometry of the current position.
    /// Since the geometry is the identity, it is here the current basis
    virtual gsBasis<T>& currentDomainBasis(gsBasis<T>& basis)
    {  return this->currentBasis(basis);  }
    virtual gsBasis<T>& currentDomainBasis(gsMultiBasis<T>& basis)
    {  return this->currentBasis(basis);  }

    /// Returns the geometry of the patch of geom in which
    /// is included the image by the geometry of the current position.
    /// Since the geometry is the identity, it is here the
    /// geometry at the current patch.
    inline gsGeometry<T>& currentDomainGeom(gsMultiBasis<T>& basis,
                                            gsFunctionSet<T>& geom)
    {
        return currentDomainGeom(getGeometry(basis, geom));
    }
    inline gsGeometry<T>& currentDomainGeom(gsBasis<T>& basis,
                                            gsFunctionSet<T>& geom)
    {
        return currentDomainGeom(getGeometry(basis, geom));
    }
    virtual gsGeometry<T>& currentDomainGeom(gsGeometry<T>& geom)
    {  return this->currentGeom(geom);  }
    virtual gsGeometry<T>& currentDomainGeom(gsMultiPatch<T>& geom)
    {  return this->currentGeom(geom);  }

}; // class gsGenGeomIterator


/*
  @brief Class gsGenGeomIteratorSimplePatch: iterator on a
  single patch geometry (see class gsGenGeomIterator).
*/
template<class T>
class gsGenGeomIteratorSimplePatch : public gsGenGeomIterator<T>
{
    typedef gsGenGeomIterator<T> Base;
public:
    /// constructor
    gsGenGeomIteratorSimplePatch(gsGeometry<T>& geom) :
    Base::gsGenGeomIterator(geom.basis()), m_geom(geom)
    {  }

    const gsVector<T> imageUpperCorner() const
    {
        return m_geom.eval(this->upperCorner());
    }
    const gsVector<T> imageLowerCorner() const
    {
        return m_geom.eval(this->lowerCorner());
    }
    const gsVector<T> imageCenterPoint() const
    {
        return m_geom.eval(this->centerPoint());
    }

    const void eval_into(gsMatrix<T>& points,
                         gsMatrix<T>& res) const
    {
        m_geom.eval_into(points, res);
    }

    /*  /// res: cols() = number of points giving the directions (d+1)
    ///      rows() = dimension of the ambient space
    void imageElementaryCell_into(gsMatrix<T> &res) const
    {
        unsigned d = m_geom.parDim();
        res.col(0) = imageLowerCorner();
        gsVector<T> diag = this->upperCorner() - this->lowerCorner();

        for(unsigned i = 0;i < d;i++)
        {
            gsVector<T> dir = this->lowerCorner();
            dir.row(i) += diag.row(i);
            res.col(i + 1) = m_geom.eval(dir);
        }
        }*/

    const unsigned basisDomDim() const
    {  return m_geom.targetDim();  }


    gsBasis<T>& currentDomainBasis(gsBasis<T>& basis)
    {  return basis;  }

    /// The domain can only have one patch in the current implementation.
    gsBasis<T>& currentDomainBasis(gsMultiBasis<T>& basis)
    {  return basis[0];  }


    gsGeometry<T>& currentDomainGeom(gsGeometry<T>& geom)
    {  return geom;  }

    /// The domain can only have one patch in the current implementation.
    gsGeometry<T>& currentDomainGeom(gsMultiBasis<T>& geom)
    {  return geom[0];  }

protected:
    /// A geometry whose domain is the same as the domain on which
    /// is defined the iterator.
    gsGeometry<T>& m_geom;

    using Base::m_patch;
}; // class gsGenGeomIteratorSimplePatch


/*
  @brief Class gsGenGeomIteratorMultiPatch: iterator on a
  multipatch geometry (see class gsGenGeomIterator).
*/
template<class T>
class gsGenGeomIteratorMultiPatch : public gsGenGeomIterator<T>
{
    typedef gsGenDomainIterator<T> Base;

public:
    /// constructor
    gsGenGeomIteratorMultiPatch(gsMultiPatch<T>& geom) :
    gsGenGeomIterator<T>::gsGenGeomIterator(geom), m_geom(geom)
    {  }

    const gsVector<T> imageUpperCorner() const
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_geom[m_patch].eval(this->upperCorner());
    }
    const gsVector<T> imageLowerCorner() const
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_geom[m_patch].eval(this->lowerCorner());
    }
    const gsVector<T> imageCenterPoint() const
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        return m_geom[m_patch].eval(this->centerPoint());
    }

    const void eval_into(gsMatrix<T>& points,
                         gsMatrix<T>& res) const
    {
        GISMO_ASSERT(m_is_good, "Error: the iteration is over");
        m_geom[m_patch].eval_into(points, res);
    }

    const unsigned basisDomDim() const
    {  return m_geom[0].targetDim();  }

    gsBasis<T>& currentDomainBasis(gsBasis<T>& basis)
    {  return basis;  }

    /// The domain can only have one patch in the current implementation.
    gsBasis<T>& currentDomainBasis(gsMultiBasis<T>& basis)
    {  return basis[0];  }


    gsGeometry<T>& currentDomainBasis(gsGeometry<T>& geom)
    {  return geom;  }

    /// The domain can only have one patch in the current implementation.
    gsGeometry<T>& currentDomainBasis(gsMultiPatch<T>& geom)
    {  return geom[0];  }

protected:
    using Base::m_patch;
    using Base::m_is_good;

    /// A multipatch geometry whose domain is the same as the
    /// domain on which is defined the iterator.
    gsMultiPatch<T>& m_geom;
}; // class gsGenGeomIteratorMultiPatch


}

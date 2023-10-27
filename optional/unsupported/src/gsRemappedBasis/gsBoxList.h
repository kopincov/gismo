/** @file gsBoxList.h

    @brief Implementation of gsBoxList, which is used in gsSelector.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsRemappedBasis/gsVectorUtils.h>

namespace gismo
{


// some utils for boxes

template <typename MatrixT>
inline bool notEmpty (const MatrixT &box)
{
    GISMO_ASSERT(box.cols()==2,"Boxes are represented by a matrix with 2 columns");
    return (box.col(0).array()<box.col(1).array()).all();
}

template <typename MatrixT>
inline bool empty (const MatrixT &box)
{
    GISMO_ASSERT(box.cols()==2,"Boxes are represented by a matrix with 2 columns");
    return !notEmpty(box);
}

template <typename MatrixT, typename MatrixT2>
inline bool intersect (const MatrixT &box1, const MatrixT2 &box2)
{
    GISMO_ASSERT(box1.cols()==2 && box2.cols()==2,"Boxes are represented by a matrix with 2 columns");
    GISMO_ASSERT(box1.rows()==box2.rows(),     "Dimension mismatch");
    for (index_t r=0; r<box1.rows(); ++r)
        if (box1(r,0)>=box2(r,1) || box2(r,0)>=box1(r,1))
            return false;
    return true;
}

template <typename MatrixT, typename MatrixT2>
inline bool intersectClosed (const MatrixT &box1, const MatrixT2 &box2)
{
    GISMO_ASSERT(box1.cols()==2 && box2.cols()==2,"Boxes are represented by a matrix with 2 columns");
    GISMO_ASSERT(box1.rows()==box2.rows(),     "Dimension mismatch");
    for (index_t r=0; r<box1.rows(); ++r)
        if (box1(r,0)>box2(r,1) || box2(r,0)>box1(r,1))
            return false;
    return true;
}

// return true if box1 is contained in box2
template <typename MatrixT, typename MatrixT2>
inline bool contained (const MatrixT &box1, const MatrixT2 &box2)
{
    GISMO_ASSERT(box1.cols()==2 && box2.cols()==2,"Boxes are represented by a matrix with 2 columns");
    GISMO_ASSERT(box1.rows()==box2.rows(),     "Dimension mismatch");
    return empty(box1) || ( (box1.col(0).array()>=box2.col(0).array() && box1.col(1).array()<=box2.col(1).array()).all());
}

/**
* @brief Checks whether @a point is inside @a box.
* Note that the boxes are closed from left and open from right.
* @param box is given as a matrix with d rows and 2 columns
* @param point is given as a matrix with d rows and 1 column
* @return true iff @a point is inside @a box
*/
template <typename MatrixT, typename MatrixT2>
static inline bool contains ( const MatrixT  &box,
                              const MatrixT2 &point)
{
    GISMO_ASSERT(box.cols()==2,"Boxes are represented by a matrix with 2 columns");
    GISMO_ASSERT(point.cols()==1,"Points are represented by a matrix with 1 columns");
    GISMO_ASSERT(box.rows()==point.rows(),     "Dimension mismatch");
    return (box.col(0).array()<=point.col(0).array()).all()
        && (box.col(1).array()> point.col(0).array()).all();
}

/** A list of boxes that together form a partition of the domain of
 * interest. For each box additional information may be stored. It is
 * used in gsSelector.*/
struct GISMO_EXPORT gsBoxList
{
    typedef unsigned basisIdT;
    typedef int directionT;
    typedef memory::shared_ptr<gsFunctionSet<real_t> > basisPtr;
    
public: // Constructors
    gsBoxList(directionT domainDim, bool basisId=true)
        : m_dim(domainDim), m_haveBasisId(basisId)
    {}

    /// Converts all the boxes in \a boxes (which are in gsHTensorBasis::refine() format) to the format
    /// used here and adds them to the box list.
    gsBoxList(const std::vector<basisPtr> &bases,
               directionT domainDim,
               const std::vector<index_t> &boxes );
    template <typename baisisT>
    gsBoxList( const std::vector<baisisT*> &bases,
               directionT domainDim,
               const std::vector<index_t> &boxes );
public: // Getters
    /// Returns the number of boxes in the boxList.
    size_t                  size    (        ) const {return m_boxes.size()/2/m_dim;}

    /// Returns the \a p-th box as a matrix.
    gsAsMatrix<real_t>      box     (size_t p)       {return gsAsMatrix<>(&m_boxes[m_dim*2*p],m_dim,2);}

    /// Returns the \a p-th box as a constant matrix.
    gsAsConstMatrix<real_t> box     (size_t p) const {return gsAsConstMatrix<>(&m_boxes[m_dim*2*p],m_dim,2);}

    /// Returns the basisId (i.e., the level, in case of a nested hierarchy) of the \a p-th box.
    basisIdT                basisId (size_t p) const {return m_basisId[p];}

    /// Returns a reference to basisId (i.e., the level, in case of a nested hierarchy) of the \a p-th box.
    basisIdT&               basisId (size_t p)       {return m_basisId[p];}

    /// Returns the minimum of the basisIds in the boxList.
    basisIdT                minId   (        ) const {return size()>0? *std::min_element(m_basisId.begin(),m_basisId.end()):0;}

    /// Returns the maximum of the basisIds in the boxList.
    basisIdT                maxId   (        ) const {return size()>0? *std::max_element(m_basisId.begin(),m_basisId.end()):0;}

    /// Returns the maximum of the basisIds of the boxes that intersect bbox
    basisIdT                maxId   (const gsMatrix<real_t> &bbox) const;

    /// Returns the dimension of the domain.
    int                     domainDim     (        ) const {return m_dim;}

public: // Modifiers

    /// Appends box given by \a box of "level" \a basisId.
    void append(const gsMatrix<real_t> &box, basisIdT basisId );

    /// Appends the given boxes to the box list.
    void append(const gsBoxList &boxes );

    /// remove a box specified by its index from the list
    void                   remove (size_t index);

    /// Do a set difference
    gsBoxList & operator-=(const gsMatrix<real_t> &other);

    inline gsBoxList difference( const gsMatrix<real_t> &box, const gsMatrix<real_t> &toRemove, basisIdT id);

    /// Do a set difference, the levels of other are ignored !!
    gsBoxList & operator-=(const gsBoxList &other);

    /**
     * @brief Returns the basisId of the basis capable of depicting all the basis functions
     * (from various basisId) active at the @a point.
     * @param point is a column of coordinates.
     * @return basisId associated to @a point.
    */
    std::vector<basisIdT> getBasisAt(const gsMatrix<real_t> point ) const;

    /**
     * @brief Checks whether @a m_boxes is a valid list.
     * @return always true; the check is through GISMO_ASSERT.
     */
    bool check() const;

    /*** Converts all the boxes in the boxList to the format used in gsHTensorBasis::refineElements.
    * Id est, each box is given by 2*dim + 1 items in an std::vector<unsigned>, where the first item
    * is the level of the box, followed by the lower left and upper right corners,
    * each given by the knot indices of the level of that particular box.
    * \param hTensorBasis THB basis so that we know the knot vectors and can determine the indices;
    * \param[out] result is cleared and the boxes shall be written here.
    */
    void toRefineElementFormat( const std::vector<basisPtr> &bases,
                                std::vector<index_t> &result ) const;
    template <typename basisT>
    void toRefineElementFormat( const std::vector<basisT*> &bases,
                                std::vector<index_t> &result ) const;



    /** Converts one box to the G+SMo format (see asGismoBoxAll).
    * \param hTensorBasis THB basis because of the indices;
    * \param p id of the box in the list;
    * \param[out] boxes : the results are appended (!) here.
    */
    static void toRefineElementFormat(const gsFunctionSet<real_t> *bases,
                                      basisIdT            lvl,
                                      gsAsConstMatrix<real_t>  box,
                                      std::vector<index_t> &boxes);


    /** Converts the box, given as a subset of a std::vector<unsigned>
     * given by \a boxBegin in the gsHTensorBasis::refine()
     * format to the gsBoxList format and write it to \a result
     * and its level to \a level. */
    void boxFromRefineElementFormat(const std::vector<basisPtr> &bases,
                                     directionT domainDim,
                                     std::vector<index_t>::const_iterator inputBegin,
                                     gsMatrix<real_t> &result,
                                     index_t &level ) const;

protected: // Members
    /**
        Stores the coordinates of the boxes. The coordinates of each box
        are stored in the same format used by gsAsMatrix and should always
        be accessed that way by using the function box(index).
    **/
    std::vector<real_t> m_boxes;

    /**
        Store sthe dimension of the stored boxes. This is used in order to
        know the total number of coordinates per box (2 m_dim).
    **/
    directionT          m_dim;

    /**
        Do we have coefficients or not?
    **/
    bool                m_haveBasisId;

    /// BasisId's (i.e., levels in case of a nested hierarchy) of the boxes.
    std::vector<basisIdT> m_basisId;
};

GISMO_EXPORT std::ostream& operator<<  (std::ostream &out, const gsBoxList &boxes);

template <typename basisT>
gsBoxList::gsBoxList( const std::vector<basisT*> &bases,
                      directionT domainDim,
                      const std::vector<index_t> &boxes )
    : m_dim(domainDim), m_haveBasisId(true)
{
    const size_t blockSize=2 * domainDim + 1;
    GISMO_ASSERT( boxes.size() % blockSize == 0, "We need one or more full boxes." );

    std::vector<index_t>::const_iterator curBox=boxes.begin();
    std::vector<index_t>::const_iterator endBox=boxes.end();

    gsMatrix<real_t> box(domainDim,2);
    basisIdT         lvl;

    for( ; curBox!=endBox; curBox += blockSize )
    {
        lvl = static_cast<basisIdT>(*curBox);
        for( directionT dir = 0; dir < domainDim; ++dir )
        {
            std::vector<real_t> uKnots=bases[lvl]->knots(dir).unique();
            box( dir, 0 ) = uKnots[ curBox[dir+1]     ];
            box( dir, 1 ) = uKnots[ curBox[dir+domainDim+1] ];
        }
        append(box, lvl);
    }
}

template <typename basisT>
void gsBoxList::toRefineElementFormat( const std::vector<basisT*> &bases,
                            std::vector<index_t> &result ) const
{
    result.clear();
    for( size_t p = 0; p < size(); ++p )
        toRefineElementFormat( bases[basisId(p)], basisId(p),box(p), result );
}


} // namespace gismo

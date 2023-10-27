/** @file gsDomainMap.h

    @brief Implementation of gsDomainMap, which is used in gsSelector.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, D. Mokris
*/

#pragma once

#include <gsRemappedBasis/tree.h>
#include <gsRemappedBasis/gsBoxList.h>

namespace gismo
{


struct NodeData
{
    typedef unsigned basisIdT;
    typedef int directionT;
    directionT dir;
    real_t     par;
    basisIdT   space;
};


template <typename T>
struct void_type {
    typedef void type;
};

template <typename Type, typename _ = void>
struct selectBoxT
{
    typedef gsMatrix<real_t,-1,-1> BoxT;
};

template <typename Type>
struct selectBoxT<Type,typename void_type<typename Type::BoxT>::type>
{
    typedef typename Type::BoxT BoxT;
};


/**
   @brief The gsDomainMap class

   A domain map represents a decomposition of a bonding box in
   sub-boxes.

   For each sub-box a basisIdT is stored. This mechanism allows
   for two different usages:

   -gsRemappedBasis uses it to decompose the domain in sub-domains
    and associate a local representation of the base to each sub-domain
    in a manifold-like manner.

   -it can be used to represent arbitrary axis-aligned polytopes by
    storing 1 in the interior and 0 in the exterior. While doing so
    the structures allows for reasonably fast boolean operations.
**/
class GISMO_EXPORT gsDomainMap : public BinaryTree<NodeData>
{
public:
    typedef unsigned basisIdT;
    typedef int directionT;
    static const basisIdT invalidBasisId = basisIdT(UINT_MAX);

protected:
    gsMatrix<real_t> m_boundingBox;
public:
    gsDomainMap()
    {
        BinaryTree<NodeData>::operator [](m_root).data.space=invalidBasisId;
    }

    gsDomainMap(
        const gsMatrix<> &boundingBox
    )
        : m_boundingBox(boundingBox)
    {
        BinaryTree<NodeData>::operator [](m_root).data.space=invalidBasisId;
    }


    gsDomainMap(
            const gsBoxList  &boxes,
            const gsMatrix<> &boundingBox
            )
    {
        initFromBoxesMax(boxes, boundingBox);
    }

    gsDomainMap(
        const gsMatrix<> &boundingBox,
        const gsMatrix<> &takenBox)
    {
        gsBoxList boxes(boundingBox.rows());
        boxes.append(takenBox,1);
        this->initFromBoxesMax(boxes,boundingBox);
    }


    /**
     * @brief Returns the bounding box of the boxes in the subtree starting in
     * @param node
     */
    gsMatrix<real_t> getBoundingBox(NodeId node, const gsMatrix<real_t> &domain) const;


    gsMatrix<real_t> getBoundingBox(NodeId node) const { return getBoundingBoxOptimized<gsMatrix<real_t> >(node); }
    template <typename BoxT>
    BoxT getBoundingBoxOptimized (NodeId node) const
    {
        BoxT result=m_boundingBox;
        NodeId lst=node;
        while ( (*this)[node].parent != no_node )
        {
            node=(*this)[node].parent;
            const directionT dir=(*this)[node].data.dir;
            if (lst==(*this)[node].first)
                result(dir,1)=math::min((*this)[node].data.par,result(dir,1));
            else
                result(dir,0)=math::max((*this)[node].data.par,result(dir,0));
            lst=node;
        }
        return result;
    }
public:
    /**
     * Returns the basisId of the basis capable of depicting all the basis functions
     * active in the point @a p.
     */
    template <typename PointType>
    basisIdT getBasisAt (const PointType &p) const
    {
        NodeId curNode=getRoot();
        while ((*this)[curNode].isFork())
        {
            const NodeData &data=(*this)[curNode];
            if( data.par>p(data.dir,0) )
                curNode=(*this)[curNode].first;
            else
                curNode=(*this)[curNode].second;
        }
        const NodeData &data=(*this)[curNode];
        return data.space;
    }

    /**
     * Returns the basisId of the basis capable of depicting all the basis functions
     * active in the given points.
     * @param points , each given as a column of the matrix;
     * @param out returned basisId's (one per point).
     */
    void getBasisAt (const gsMatrix<real_t> &points, std::vector<basisIdT> &out) const;

    /**
       @brief getBoundingBox
       @return the bounding box of the domain
     */
    const gsMatrix<real_t> & getBoundingBox() const {return m_boundingBox;}
public:

    /**
     * Constructs the domainMap from a list of boxes and a bounding box.
     * Every time two boxes overlap the maximum basisId is taken on the intersection.
     * Area that is not covered by the boxes is set to have basisId 0.
     * @param boxes list of boxes
     * @param boundingBox bounding box
     */
    void initFromBoxesMax(
            const gsBoxList  &boxes,
            const gsMatrix<> &boundingBox
            );

    /**
     * Change the domainMap by setting the basisId in each point to the maximum
     * between what is assigned and the basisId-s of the containing boxes in the
     * provided box list.
     * @param boxes list of boxes
     */
    void addBoxesMax(
            const gsBoxList  &boxes
            );


    /**
     * Constructs the domainMap from a list of boxes and a bounding box.
     * Every time two boxes overlap it is an error.
     *
     * @param boxes list of boxes
     * @param boundingBox bounding box
     */
    void initFromBoxes(
            const gsBoxList  &boxes,
            const gsMatrix<> &boundingBox
            );
    /**
     * Change the domainMap by setting the basisId in each point to the minimum
     * between what is assigned and the basisId-s of the containing boxes in the
     * provided box list.
     * @param boxes list of boxes
     */
    void addBoxesMin(
            const gsBoxList  &boxes
            );


    /**
     * Constructs the domainMap from a list of boxes and a bounding box.
     * Every time two boxes overlap the maximum basisId is taken on the intersection.
     * Area that is not covered by the boxes is set to have basisId 0.
     * @param boxes list of boxes
     * @param boundingBox bounding box
     */
    void initFromBoxesMin(
            const gsBoxList  &boxes,
            const gsMatrix<> &boundingBox
            );


    /**
     * Set the full domain to the b.
     * @param boxes list of boxes
     * @param boundingBox bounding box
     */
    void initConstant (
        basisIdT b,
        const gsMatrix<> &boundingBox
        );

    /**
     * @brief Returns the domains in the tree restricted to @a boundingBox
     * as a list of boxes.
     */
    gsBoxList asBoxList(const gsMatrix<> &boundingBox) const;
    gsBoxList asBoxList() const {return asBoxList(m_boundingBox);}


    /**
      Apply the action to the subtree starting in node and intersecting
      the given boundingBox bBox.

      \param[in] node, the starting node (usually root() )
      \param[in] bBox, a matrix containing the a box defining the area on
                 which the action should be executed
      \param[in] action, an instance of a struct providing two methods with
                 the following signature:

                    void enter (gsDomainMap *, gsMatrix<real_t,dim1,dim2>, NodeId)
                    void exit  (gsDomainMap *, gsMatrix<real_t,dim1,dim2>, NodeId)

      Traverse the tree nodes that intersect whose bounding box intersects
      \a bBox and that are descendent of \a node. For each node \a cur the
      provided function action.enter is called before traversing the subtree
      of \a cur and action.exit is called after traversing the same subtree.
     */

    template <typename Action>
    bool applyToSubtree(NodeId node, const typename selectBoxT<Action>::BoxT &bBox, Action &action);
    template <typename Action>
    bool applyToSubtree(NodeId node, const typename selectBoxT<Action>::BoxT &bBox, Action &action) const;

    template <typename Action>
    bool apply(const typename selectBoxT<Action>::BoxT &bBox, Action &action);
    template <typename Action>
    bool apply(const typename selectBoxT<Action>::BoxT &bBox, Action &action) const;

    template <typename Action>
    bool apply(Action &action);
    template <typename Action>
    bool apply(Action &action)const;

public:
    // actions on the tree used with applyToSubtree

    // base

    template <int dim1=-1, int dim2=-1>
    struct actionBase;

    // get list of boxes
    template <int dim1=-1,int dim2=-1>
    struct getBoxes;

    // modification

    template <int dim1=-1, int dim2=-1>
    struct optimizeOnExitAction;

    template <int dim1=-1,int dim2=-1>
    struct setBasis;

    template <int dim1=-1,int dim2=-1>
    struct shiftIds;

    // used by refinement and unions
    template <int dim1=-1,int dim2=-1>
    struct setBasisMax;

    // used by derefinement and intersections
    template <int dim1=-1,int dim2=-1>
    struct setBasisMin;

    // used by difference
    template <int dim1=-1,int dim2=-1>
    struct setDifference;
};


GISMO_EXPORT gsDomainMap shiftIds ( const gsDomainMap & orig, int amount);

GISMO_EXPORT gsDomainMap polytopeUnion( const gsDomainMap & left,
                                        const gsDomainMap & right );

GISMO_EXPORT gsDomainMap polytopeIntersect( const gsDomainMap & left,
                                            const gsDomainMap & right );

GISMO_EXPORT gsDomainMap polytopeDifference( const gsDomainMap & left,
                                             const gsDomainMap & right );


GISMO_EXPORT  std::ostream& operator<<  (std::ostream &out, const NodeData &node);
} // namespace gismo

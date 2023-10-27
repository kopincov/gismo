

#pragma once

#include <gsRemappedBasis/gsDomainMap.h>


namespace gismo {

// NOTE THE CODE DUPLICATION FOR CONST AND NOT CONST VERSION
// NO IDEA HOW TO REMOVE IT

template <typename Action>
bool gsDomainMap::applyToSubtree(NodeId node, const typename selectBoxT<Action>::BoxT &bBox, Action &action)
{
    bool cont=action.enter( this, bBox, node);

    if ( cont && (*this)[node].isFork() )
    {
        typename selectBoxT<Action>::BoxT boxCopy=bBox;

        const directionT dir=(*this)[node].data.dir;
        const real_t     par=(*this)[node].data.par;
        const real_t     end=bBox(dir,1);

        if ( cont && (*this)[node].first != no_node && bBox(dir,0)<par )
        {
            boxCopy(dir,1)=math::min(end,par);
            cont=applyToSubtree((*this)[node].first, boxCopy , action);
            boxCopy(dir,1)=end;
        }
        if ( cont && (*this)[node].second != no_node && bBox(dir,1)>par )
        {
            boxCopy(dir,0)=math::max(bBox(dir,0),par);
            cont=applyToSubtree((*this)[node].second, boxCopy , action);
        }
    }
    return cont && action.exit( this, bBox, node);
}

// const version
template <typename Action>
bool gsDomainMap::applyToSubtree(NodeId node, const typename selectBoxT<Action>::BoxT &bBox, Action &action) const
{
    bool cont=action.enter( this, bBox, node);

    if ( cont && (*this)[node].isFork() )
    {
        typename selectBoxT<Action>::BoxT boxCopy=bBox;

        const directionT dir=(*this)[node].data.dir;
        const real_t     par=(*this)[node].data.par;
        const real_t     end=bBox(dir,1);

        if ( cont && (*this)[node].first != no_node && bBox(dir,0)<par )
        {
            boxCopy(dir,1)=math::min(end,par);
            cont=applyToSubtree((*this)[node].first, boxCopy , action);
            boxCopy(dir,1)=end;
        }
        if ( cont && (*this)[node].second != no_node && bBox(dir,1)>par )
        {
            boxCopy(dir,0)=math::max(bBox(dir,0),par);
            cont=applyToSubtree((*this)[node].second, boxCopy , action);
        }
    }
    return cont && action.exit( this, bBox, node);
}

template <typename Action>
bool gsDomainMap::apply(const typename selectBoxT<Action>::BoxT &bBox, Action &action) const { return applyToSubtree(m_root,bBox,action); }
template <typename Action>
bool gsDomainMap::apply(const typename selectBoxT<Action>::BoxT &bBox, Action &action)       { return applyToSubtree(m_root,bBox,action); }


template <typename Action>
bool gsDomainMap::apply(Action &action)       { return applyToSubtree(m_root,m_boundingBox,action); }
template <typename Action>
bool gsDomainMap::apply(Action &action)const  { return applyToSubtree(m_root,m_boundingBox,action); }


// Actions

template <int dim1, int dim2>
struct gsDomainMap::actionBase
{
    typedef gsMatrix<real_t,dim1,dim2> BoxT;
    bool enter (const gsDomainMap *, const BoxT & , NodeId ) {return true;}
    bool exit  (const gsDomainMap *, const BoxT & , NodeId ) {return true;}

    static void split (gsDomainMap *map,
                       NodeId node,
                       const BoxT& box,
                       const BoxT& leafBox)
    {
        const directionT domainDim = static_cast<directionT>(box.rows());

        NodeId first =map->append(node,true);
        NodeId second=map->append(node,false);
        (*map)[first].data.space  =(*map)[node].data.space;
        (*map)[second].data.space =(*map)[node].data.space;

        for (directionT dir=0; dir<domainDim;++dir)
        {
            if (leafBox(dir,0)<box(dir,0))
            {
                (*map)[node].data.dir=dir;
                (*map)[node].data.par=box(dir,0);
                break;
            }
            else if (leafBox(dir,1)>box(dir,1))
            {
                (*map)[node].data.dir=dir;
                (*map)[node].data.par=box(dir,1);
                break;
            }
        }
    }
};

template <int dim1, int dim2>
struct gsDomainMap::optimizeOnExitAction : public actionBase<dim1,dim2>
{
    typedef gsMatrix<real_t,dim1,dim2> BoxT;
    /// Split in the first useful direction.

    bool exit (gsDomainMap *map,
               const BoxT & /*box*/,
               NodeId       node )
    {
        if ((*map)[node].isLeaf())
            return true;
        // compact the tree by removing unnecessary splits and leaves
        // in the subtree changed
        removeCommonDirSplit2(map, node);
        return true;

        //        if ( removeConstSplit(map, node) )
        //            return true;
        //        while ( removeInterSplit(map, node, &Node::first) )
        //            ;
        //        while ( removeInterSplit(map, node, &Node::second) )
        //            ;
        //        return true;
    }


    static inline bool compareSubtree (gsDomainMap &map, NodeId node1, NodeId node2, int maxNodes)
    {
        if ( map[node1].isFork() && map[node2].isFork() && maxNodes>0)
        {
            if (    map[node1].data.par != map[node2].data.par
                    ||  map[node1].data.dir != map[node2].data.dir
                    ||  ( (map[node1].first ==no_node) != (map[node2].first ==no_node ) )
                    ||  ( (map[node1].second==no_node) != (map[node2].second==no_node ) )
                    )
                return false;
            if ( ( map[node1].first==no_node) || compareSubtree (map, map[node1].first, map[node2].first, maxNodes-1) )
                return ( ( map[node1].second==no_node) || compareSubtree (map, map[node1].second, map[node2].second, maxNodes-1) );
            else return false;
        }
        else
            return map[node1].isLeaf()
                && map[node2].isLeaf()
                && (map[node1].data.space == map[node2].data.space);
    }



    inline static bool removeCommonDirSplit2 (gsDomainMap *map, NodeId node)
    {
        NodeId lastB=(*map)[node].first;  // last before
        NodeId frstA=(*map)[node].second; // first after
        NodeId splitB=node;
        NodeId splitA=node;
        bool changed=false;

        if (lastB==no_node ||frstA==no_node)
            return changed;
        const directionT dir=(*map)[node].data.dir;

        while ( (*map)[lastB].second != no_node && (*map)[lastB].data.dir==dir)
            lastB=(*map)[lastB].second;
        while ( (*map)[frstA].first  != no_node && (*map)[frstA].data.dir==dir)
            frstA=(*map)[frstA].first;
        while ( compareSubtree (*map, frstA, lastB, 10) )
        {
            changed = true;
            splitB = (*map)[lastB].parent;
            splitA = (*map)[frstA].parent;
            if (splitB==splitA) // constant split on node
            {
                map->remove(lastB);
                map->detach(frstA);
                map->swap(frstA,splitA);
                break; // can not simplify further
            }
            else if ( splitB != node  )
            {
                (*map)[node].data.par = (*map)[splitB].data.par;
                lastB=(*map)[splitB].first;
                map->swap  (splitB, lastB);
                map->remove(splitB);
                while ( (*map)[lastB].second != no_node && (*map)[lastB].data.dir==dir )
                    lastB=(*map)[lastB].second;
            }
            else
            {
                (*map)[node].data.par = (*map)[splitA].data.par;
                frstA=(*map)[splitA].second;
                map->swap(splitA, frstA);
                map->remove(splitA);
                while ( (*map)[frstA].first  != no_node && (*map)[frstA].data.dir==dir)
                    frstA=(*map)[frstA].first;
            }
        }
        return changed;
    }


    inline static bool removeCommonDirSplit (gsDomainMap *map, NodeId node)
    {
        NodeId lastB=(*map)[node].first;  // last before
        NodeId frstA=(*map)[node].second; // first after
        NodeId splitB=node;
        NodeId splitA=node;
        bool changed=false;

        if (lastB==no_node ||frstA==no_node)
            return changed;
        const directionT dir=(*map)[node].data.dir;

        while ( (*map)[lastB].second != no_node && (*map)[lastB].data.dir==dir)
            lastB=(*map)[lastB].second;
        if (!(*map)[lastB].isLeaf())
            return changed;
        while ( (*map)[frstA].first  != no_node && (*map)[frstA].data.dir==dir)
            frstA=(*map)[frstA].first;
        if (!(*map)[frstA].isLeaf())
            return changed;
        while ( (*map)[frstA].data.space == (*map)[lastB].data.space )
        {
            changed = true;
            splitB = (*map)[lastB].parent;
            splitA = (*map)[frstA].parent;
            if (splitB==splitA) // constant split on node
            {
                (*map)[node].data.space=(*map)[lastB].data.space;
                map->remove(lastB);
                map->remove(frstA);
                break; // can not simplify further
            }
            else if ( splitB != node  )
            {
                (*map)[node].data.par = (*map)[splitB].data.par;
                lastB=(*map)[splitB].first;
                map->swap  (splitB, lastB);
                map->remove(splitB);
                while ( (*map)[lastB].second != no_node && (*map)[lastB].data.dir==dir )
                    lastB=(*map)[lastB].second;
                if (!(*map)[lastB].isLeaf())
                    return changed;
            }
            else
            {
                (*map)[node].data.par = (*map)[splitA].data.par;
                frstA=(*map)[splitA].second;
                map->swap(splitA, frstA);
                map->remove(splitA);
                while ( (*map)[frstA].first  != no_node && (*map)[frstA].data.dir==dir)
                    frstA=(*map)[frstA].first;
                if (!(*map)[frstA].isLeaf())
                    return changed;
            }
        }
        return changed;
    }

    inline static bool removeInterSplit (gsDomainMap *map, NodeId node, NodeId Node::*leafDir)
    {
        NodeId Node::*forkDir =  leafDir==&Node::first ? &Node::second : &Node::first;

        const NodeId  leaf = (*map)[node].*leafDir;
        const NodeId  fork = (*map)[node].*forkDir;
        if (leaf==no_node || fork==no_node ||  (*map)[leaf].isFork())
            return false;
        const directionT dir =  (*map)[node].data.dir;

        NodeId nephew=fork;
        NodeId nephewPar=node;
        while ( nephew != no_node && (*map)[nephew].isFork() )
        {
            if ( (*map)[nephew].data.dir != dir ) // split in different direction -> abort
                return false;
            else
            {
                nephewPar = nephew;
                nephew    = (*map)[nephew].*leafDir;
            }
        }
        if (nephew == no_node || (*map)[nephew].data.space!=(*map)[leaf].data.space )
            return false;

        (*map)[node].data.par=(*map)[nephewPar].data.par;
        map->swap(nephewPar,(*map)[nephewPar].*forkDir);
        map->remove(nephewPar);
        return true;
    }

    inline static bool removeConstSplit (gsDomainMap *map, NodeId node)
    {
        const NodeId  par = (*map)[node].parent;
        if (par==no_node)
            return false;
        const NodeId  fst = (*map)[par].first;
        const NodeId  snd = (*map)[par].second;
        if (fst==no_node || snd==no_node ||(*map)[fst].isFork() || (*map)[snd].isFork()|| (*map)[fst].data.space!=(*map)[snd].data.space )
            return false;
        (*map)[node].data.space=(*map)[fst].data.space;
        map->remove(fst);
        map->remove(snd);
        return true;
    }
};

template <int dim1, int dim2>
struct gsDomainMap::setBasis : public gsDomainMap::optimizeOnExitAction<dim1,dim2>
{
    typedef gsMatrix<real_t,dim1,dim2> BoxT;
    basisIdT m_bas;

    setBasis(basisIdT basis)
        : m_bas(basis)
    {}
    bool enter (gsDomainMap *map, const BoxT &box, NodeId node )
    {
        if ((*map)[node].isFork() )
            return true;
        BoxT leafBox=map->getBoundingBoxOptimized<BoxT>(node);
        if (contained(leafBox, box))
            (*map)[node].data.space=m_bas;
        else
            actionBase<dim1,dim2>::split(map, node, box, leafBox);
        return true;
    }

};

template <int dim1, int dim2>
struct gsDomainMap::setBasisMax : public gsDomainMap::setBasis<dim1,dim2>
{
    typedef gsMatrix<real_t,dim1,dim2> BoxT;
    using gsDomainMap::setBasis<dim1,dim2>::m_bas;

    setBasisMax(basisIdT basis)
        : setBasis<dim1,dim2>(basis)
    {}
    bool enter (gsDomainMap *map, const BoxT &box, NodeId node )
    {
        // either increase the level or split depending if the
        // current leaf is fully contained in the box or just
        // intersect
        if ((*map)[node].isFork() || m_bas<=(*map)[node].data.space)
            return true;

        BoxT leafBox=map->getBoundingBoxOptimized<BoxT>(node);
        if (contained(leafBox, box))
            (*map)[node].data.space=m_bas;
        else
            actionBase<dim1,dim2>::split(map, node, box, leafBox);
        return true;
    }
};

template <int dim1, int dim2>
struct gsDomainMap::setBasisMin : public gsDomainMap::setBasis<dim1,dim2>
{
    typedef gsMatrix<real_t,dim1,dim2> BoxT;
    using gsDomainMap::setBasis<dim1,dim2>::m_bas;
    setBasisMin(basisIdT basis)
        : setBasis<dim1,dim2>(basis)
    {}
    bool enter (gsDomainMap *map, const BoxT &box, NodeId node )
    {
        // either decrease the level or split depending if the
        // current leaf is fully contained in the box or just
        // intersect
        if ((*map)[node].isFork() || m_bas>=(*map)[node].data.space)
            return true;

        BoxT leafBox=map->getBoundingBoxOptimized<BoxT>(node);
        if (contained(leafBox, box))
            (*map)[node].data.space=m_bas;
        else
            actionBase<dim1,dim2>::split(map, node, box, leafBox);
        return true;
    }
};




template <int dim1, int dim2>
struct gsDomainMap::setDifference : public gsDomainMap::setBasis<dim1,dim2>
{
    typedef gsMatrix<real_t,dim1,dim2> BoxT;
    using gsDomainMap::setBasis<dim1,dim2>::m_bas;

    setDifference(basisIdT a) : setBasis<dim1,dim2>(a) {}

    bool enter(gsDomainMap *map,
               const BoxT &box,
               NodeId node )
    {
        if (m_bas == 0)
            return false;
        if ((*map)[node].isFork() )
            return true;

        BoxT leafBox=map->getBoundingBoxOptimized<BoxT>(node);
        if (contained(leafBox, box))
            (*map)[node].data.space = 0;
        else
            actionBase<dim1,dim2>::split(map, node, box, leafBox);
        return true;
    }
};

template <int dim1, int dim2>
struct gsDomainMap::shiftIds : public gsDomainMap::actionBase<dim1,dim2>
{
    typedef gsMatrix<real_t,dim1,dim2> BoxT;

    int m_shift;

    shiftIds(int shift)
        : m_shift(shift)
    {}
    bool enter (gsDomainMap *map,
                const BoxT & /*box*/,
                NodeId       node )
    {
        if ((*map)[node].isFork() )
            return true;
        (*map)[node].data.space+=m_shift;
        return true;
    }
};


template <template <int,int> class Action,int dim1, int dim2>
struct applyToOther : public gsDomainMap::actionBase<dim1,dim2>
{
    typedef gsMatrix<real_t,dim1,dim2> BoxT;
    gsDomainMap &target;

    applyToOther(gsDomainMap &other)
        : target(other)
    {}

    bool enter  (const gsDomainMap *map, const gsMatrix<real_t,dim1,dim2> &bbox, gsDomainMap::NodeId &node )
    {
        if ((*map)[node].isFork())
            return true;
        typename selectBoxT<Action<dim1,dim2> >::BoxT boxC=bbox;
        Action<dim1,dim2> action((*map)[node].data.space);
        target.apply(boxC, action);
        return true;
    }
};


template <int dim1, int dim2>
struct gsDomainMap::getBoxes : public actionBase<dim1,dim2>
{
    typedef gsMatrix<real_t,dim1,dim2> BoxT;
    gsBoxList &m_list;

    getBoxes(gsBoxList &list)
        : m_list(list)
    {}

    bool enter  (const gsDomainMap *map, const BoxT &bbox, NodeId &node )
    {
        if ((*map)[node].isFork())
            return true;
        m_list.append(bbox,(*map)[node].data.space);
        return true;
    }
};


} //namespace gismo

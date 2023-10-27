
/********************************************
  Graph Representation using adjacency list.

  based on cgt
*********************************************/

#pragma once


//#include <gsExternal/cgt/graph.h>
# include <ostream>

namespace gismo {

template<class NodeType, class EdgeType> //, class GraphType=cgt::_Directed>
class gsGraph // : public cgt::graph<NodeType,EdgeType,GraphType>
{
public:

    // Initialize empty graph
    gsGraph() { };

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const ;

    size_t size() const
    {
        size_t k(0);
        //typename cgt::graph<NodeType,EdgeType,GraphType>::const_iterator it;
        //for ( it =  this->begin(); it!= this->end(); ++it ) ++k;
        return k;
    };

    NodeType operator[] (size_t const & i) const
    {
        //typename cgt::graph<NodeType,EdgeType,GraphType>::const_iterator it=  this->begin();
        //for ( size_t k = 0; k<i; ++k, ++it );
        //    return  it->value();
        return NodeType();
    };
};


template<class NodeType, class EdgeType>
std::ostream & 
gsGraph<NodeType,EdgeType>::print(std::ostream &os) const 
{
    os << "gsGraph ("<<(this->is_directed () ? "directed":"undirected")<<"): " << std::endl;

    // typename gsGraph<NodeType,EdgeType,GraphType>::const_iterator itn;
    // for (itn = this->begin (); itn != this->end (); ++itn)
    // {
    //     const typename gsGraph<NodeType,EdgeType,GraphType>::node&     n = *itn;
    //     const typename gsGraph<NodeType,EdgeType,GraphType>::vertex&  v = n.vertex ();

    //     os << "vertex: " << v.value () << std::endl;

    //     const typename gsGraph<NodeType,EdgeType,GraphType>::adjlist &adjList = n.adjlist ();

    //     typename gsGraph<NodeType,EdgeType,GraphType>::adjlist::const_iterator itadj;
    //     for (itadj = adjList.begin (); itadj != adjList.end(); ++itadj)
    //     {
    //         const typename gsGraph<NodeType,EdgeType,GraphType>::edge&     e   = itadj->edge ();
    //         const typename gsGraph<NodeType,EdgeType,GraphType>::vertex&  v1  = e.v1 ();
    //         const typename gsGraph<NodeType,EdgeType,GraphType>::vertex&  v2  = e.v2 ();

    //         os << "  edge (" << e.value () << ", " << v1.value () << ", " << v2.value () << ")";
    //     }
    // }
    return os;
};


/// Print (as string) operator to be used by all derived classes
template<class NodeType, class EdgeType>
std::ostream &operator<<(std::ostream &os, const gsGraph<NodeType,EdgeType>& b)
{return b.print(os); };

} ; // namespace gismo

